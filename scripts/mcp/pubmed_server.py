#!/usr/bin/env python3
"""
PubMed MCP Server for InterSubMod Project

This server provides PubMed search capabilities via the Model Context Protocol (MCP).
It uses NCBI E-utilities API to search and retrieve biomedical literature.

Usage:
    The server will automatically load API key from:
    1. Environment variable NCBI_API_KEY
    2. .env file in the project root

    python pubmed_server.py

Environment Variables:
    NCBI_API_KEY: Your NCBI API key for higher rate limits
"""

import json
import os
import sys
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Any, Dict, List, Optional
import xml.etree.ElementTree as ET


def load_env_file():
    """Load environment variables from .env file."""
    # Try to find .env file in project root
    current_dir = Path(__file__).parent
    project_root = current_dir.parent.parent  # scripts/mcp -> scripts -> project_root
    env_file = project_root / ".env"

    if env_file.exists():
        with open(env_file, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, value = line.split("=", 1)
                    key = key.strip()
                    value = value.strip().strip('"').strip("'")
                    if key and value and key not in os.environ:
                        os.environ[key] = value


# Load .env file on import
load_env_file()


class PubMedClient:
    """Client for interacting with PubMed E-utilities API."""

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("NCBI_API_KEY", "")

    def _make_request(self, endpoint: str, params: Dict[str, str]) -> str:
        """Make a request to the E-utilities API."""
        if self.api_key:
            params["api_key"] = self.api_key

        params["retmode"] = "xml"
        query_string = urllib.parse.urlencode(params)
        url = f"{self.BASE_URL}/{endpoint}?{query_string}"

        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                return response.read().decode("utf-8")
        except Exception as e:
            return f"Error: {str(e)}"

    def search(self, query: str, max_results: int = 10) -> List[str]:
        """Search PubMed and return list of PMIDs."""
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": str(max_results),
            "sort": "relevance"
        }

        response = self._make_request("esearch.fcgi", params)

        try:
            root = ET.fromstring(response)
            pmids = [id_elem.text for id_elem in root.findall(".//Id") if id_elem.text]
            return pmids
        except ET.ParseError:
            return []

    def fetch_details(self, pmids: List[str]) -> List[Dict[str, Any]]:
        """Fetch article details for given PMIDs."""
        if not pmids:
            return []

        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "rettype": "abstract"
        }

        response = self._make_request("efetch.fcgi", params)

        articles = []
        try:
            root = ET.fromstring(response)
            for article in root.findall(".//PubmedArticle"):
                article_data = self._parse_article(article)
                if article_data:
                    articles.append(article_data)
        except ET.ParseError:
            pass

        return articles

    def _parse_article(self, article_elem: ET.Element) -> Optional[Dict[str, Any]]:
        """Parse a PubmedArticle XML element into a dictionary."""
        try:
            medline = article_elem.find(".//MedlineCitation")
            if medline is None:
                return None

            pmid_elem = medline.find("PMID")
            pmid = pmid_elem.text if pmid_elem is not None else "Unknown"

            article = medline.find("Article")
            if article is None:
                return None

            # Title
            title_elem = article.find("ArticleTitle")
            title = title_elem.text if title_elem is not None else "No title"

            # Abstract
            abstract_parts = []
            abstract_elem = article.find("Abstract")
            if abstract_elem is not None:
                for text in abstract_elem.findall("AbstractText"):
                    label = text.get("Label", "")
                    content = text.text or ""
                    if label:
                        abstract_parts.append(f"{label}: {content}")
                    else:
                        abstract_parts.append(content)
            abstract = " ".join(abstract_parts) if abstract_parts else "No abstract available"

            # Authors
            authors = []
            author_list = article.find("AuthorList")
            if author_list is not None:
                for author in author_list.findall("Author"):
                    last_name = author.find("LastName")
                    fore_name = author.find("ForeName")
                    if last_name is not None:
                        name = last_name.text
                        if fore_name is not None:
                            name = f"{fore_name.text} {name}"
                        authors.append(name)

            # Journal
            journal_elem = article.find(".//Journal/Title")
            journal = journal_elem.text if journal_elem is not None else "Unknown journal"

            # Year
            year_elem = article.find(".//PubDate/Year")
            year = year_elem.text if year_elem is not None else "Unknown year"

            # DOI
            doi = None
            for id_elem in article_elem.findall(".//ArticleId"):
                if id_elem.get("IdType") == "doi":
                    doi = id_elem.text
                    break

            return {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "authors": authors[:5],  # Limit to first 5 authors
                "journal": journal,
                "year": year,
                "doi": doi,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }
        except Exception:
            return None


def format_search_results(articles: List[Dict[str, Any]]) -> str:
    """Format search results for display."""
    if not articles:
        return "No results found."

    output = []
    for i, article in enumerate(articles, 1):
        authors_str = ", ".join(article["authors"][:3])
        if len(article["authors"]) > 3:
            authors_str += " et al."

        output.append(f"""
## {i}. {article['title']}

- **Authors**: {authors_str}
- **Journal**: {article['journal']} ({article['year']})
- **PMID**: {article['pmid']}
- **URL**: {article['url']}
{f"- **DOI**: {article['doi']}" if article['doi'] else ""}

**Abstract**: {article['abstract'][:500]}{"..." if len(article['abstract']) > 500 else ""}
""")

    return "\n---\n".join(output)


def main():
    """Main function for MCP server."""
    client = PubMedClient()

    # Read from stdin for MCP protocol
    for line in sys.stdin:
        try:
            request = json.loads(line.strip())

            method = request.get("method", "")
            params = request.get("params", {})

            if method == "search":
                query = params.get("query", "")
                max_results = params.get("max_results", 10)

                pmids = client.search(query, max_results)
                articles = client.fetch_details(pmids)
                result = format_search_results(articles)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(articles)
                    }
                }

            elif method == "fetch":
                pmids = params.get("pmids", [])
                articles = client.fetch_details(pmids)
                result = format_search_results(articles)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(articles)
                    }
                }

            else:
                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "error": {
                        "code": -32601,
                        "message": f"Unknown method: {method}"
                    }
                }

            print(json.dumps(response), flush=True)

        except json.JSONDecodeError:
            continue
        except Exception as e:
            error_response = {
                "jsonrpc": "2.0",
                "id": None,
                "error": {
                    "code": -32603,
                    "message": str(e)
                }
            }
            print(json.dumps(error_response), flush=True)


if __name__ == "__main__":
    # For testing without MCP
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        client = PubMedClient()
        query = " ".join(sys.argv[2:]) if len(sys.argv) > 2 else "methylation cancer"
        print(f"Searching PubMed for: {query}")
        pmids = client.search(query, 5)
        print(f"Found {len(pmids)} results")
        articles = client.fetch_details(pmids)
        print(format_search_results(articles))
    else:
        main()
