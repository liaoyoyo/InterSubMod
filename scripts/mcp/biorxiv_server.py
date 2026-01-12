#!/usr/bin/env python3
"""
bioRxiv MCP Server for InterSubMod Project

This server provides bioRxiv preprint search capabilities via the Model Context Protocol (MCP).
It uses the bioRxiv API to search and retrieve preprint papers.

Usage:
    python biorxiv_server.py

    Test mode:
    python biorxiv_server.py --test "methylation cancer"
"""

import json
import sys
import urllib.request
import urllib.parse
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional


class BioRxivClient:
    """Client for interacting with bioRxiv API."""

    BASE_URL = "https://api.biorxiv.org"

    def __init__(self):
        pass

    def _make_request(self, endpoint: str) -> Optional[Dict]:
        """Make a request to the bioRxiv API."""
        url = f"{self.BASE_URL}/{endpoint}"

        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "InterSubMod-MCP/1.0"}
            )
            with urllib.request.urlopen(req, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as e:
            return {"error": str(e)}

    def search_by_subject(self, subject: str, interval: str = "30d", cursor: int = 0) -> List[Dict[str, Any]]:
        """Search bioRxiv by subject area within a date interval.

        Args:
            subject: Subject area (e.g., "bioinformatics", "genomics", "cancer-biology")
            interval: Date interval (e.g., "30d" for last 30 days)
            cursor: Pagination cursor

        Returns:
            List of paper dictionaries
        """
        # Calculate date range
        end_date = datetime.now()
        days = int(interval.replace("d", ""))
        start_date = end_date - timedelta(days=days)

        start_str = start_date.strftime("%Y-%m-%d")
        end_str = end_date.strftime("%Y-%m-%d")

        # bioRxiv API endpoint for content by date
        endpoint = f"details/biorxiv/{start_str}/{end_str}/{cursor}"

        response = self._make_request(endpoint)

        if not response or "collection" not in response:
            return []

        # Filter by subject if specified
        papers = response.get("collection", [])
        if subject:
            subject_lower = subject.lower()
            papers = [
                p for p in papers
                if subject_lower in p.get("category", "").lower()
                or subject_lower in p.get("title", "").lower()
                or subject_lower in p.get("abstract", "").lower()
            ]

        return papers[:10]  # Limit to 10 results

    def search_by_keyword(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """Search bioRxiv by keyword in recent papers.

        Args:
            query: Search query
            max_results: Maximum number of results to return

        Returns:
            List of paper dictionaries
        """
        # Get recent papers (last 60 days) and filter by keyword
        end_date = datetime.now()
        start_date = end_date - timedelta(days=60)

        start_str = start_date.strftime("%Y-%m-%d")
        end_str = end_date.strftime("%Y-%m-%d")

        endpoint = f"details/biorxiv/{start_str}/{end_str}/0"

        response = self._make_request(endpoint)

        if not response or "collection" not in response:
            return []

        papers = response.get("collection", [])

        # Filter by keyword
        query_lower = query.lower()
        query_terms = query_lower.split()

        filtered_papers = []
        for paper in papers:
            title = paper.get("title", "").lower()
            abstract = paper.get("abstract", "").lower()
            category = paper.get("category", "").lower()

            # Check if any query term matches
            if any(term in title or term in abstract or term in category for term in query_terms):
                filtered_papers.append(paper)

            if len(filtered_papers) >= max_results:
                break

        return filtered_papers

    def get_paper_details(self, doi: str) -> Optional[Dict[str, Any]]:
        """Get details for a specific paper by DOI.

        Args:
            doi: Paper DOI (e.g., "10.1101/2024.01.01.000000")

        Returns:
            Paper details dictionary or None
        """
        # Extract the bioRxiv-specific part of DOI
        if "10.1101/" in doi:
            paper_id = doi.split("10.1101/")[1]
        else:
            paper_id = doi

        endpoint = f"details/biorxiv/10.1101/{paper_id}/na/na"

        response = self._make_request(endpoint)

        if response and "collection" in response and response["collection"]:
            return response["collection"][0]

        return None


def format_search_results(papers: List[Dict[str, Any]]) -> str:
    """Format search results for display."""
    if not papers:
        return "No results found."

    output = []
    for i, paper in enumerate(papers, 1):
        # Extract authors
        authors = paper.get("authors", "Unknown authors")
        if isinstance(authors, str) and len(authors) > 100:
            authors = authors[:100] + "..."

        # Extract abstract
        abstract = paper.get("abstract", "No abstract available")
        if len(abstract) > 500:
            abstract = abstract[:500] + "..."

        # Format date
        date = paper.get("date", "Unknown date")

        output.append(f"""
## {i}. {paper.get('title', 'No title')}

- **Authors**: {authors}
- **Category**: {paper.get('category', 'Unknown')}
- **Date**: {date}
- **DOI**: {paper.get('doi', 'N/A')}
- **URL**: https://www.biorxiv.org/content/{paper.get('doi', '')}

**Abstract**: {abstract}
""")

    return "\n---\n".join(output)


def main():
    """Main function for MCP server."""
    client = BioRxivClient()

    # Read from stdin for MCP protocol
    for line in sys.stdin:
        try:
            request = json.loads(line.strip())

            method = request.get("method", "")
            params = request.get("params", {})

            if method == "search":
                query = params.get("query", "")
                max_results = params.get("max_results", 10)

                papers = client.search_by_keyword(query, max_results)
                result = format_search_results(papers)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(papers)
                    }
                }

            elif method == "search_by_subject":
                subject = params.get("subject", "")
                interval = params.get("interval", "30d")

                papers = client.search_by_subject(subject, interval)
                result = format_search_results(papers)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(papers)
                    }
                }

            elif method == "get_paper":
                doi = params.get("doi", "")
                paper = client.get_paper_details(doi)

                if paper:
                    result = format_search_results([paper])
                else:
                    result = "Paper not found."

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": 1 if paper else 0
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
        client = BioRxivClient()
        query = " ".join(sys.argv[2:]) if len(sys.argv) > 2 else "methylation cancer"
        print(f"Searching bioRxiv for: {query}")
        papers = client.search_by_keyword(query, 5)
        print(f"Found {len(papers)} results")
        print(format_search_results(papers))
    else:
        main()
