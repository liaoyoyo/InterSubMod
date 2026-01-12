#!/usr/bin/env python3
"""
Ensembl MCP Server for InterSubMod Project

This server provides Ensembl gene annotation capabilities via the Model Context Protocol (MCP).
It uses the Ensembl REST API for gene queries, coordinate lookups, and variant annotations.

Usage:
    python ensembl_server.py

    Test mode:
    python ensembl_server.py --test "BRCA1"
"""

import json
import sys
import urllib.request
import urllib.parse
from typing import Any, Dict, List, Optional


class EnsemblClient:
    """Client for interacting with Ensembl REST API."""

    BASE_URL = "https://rest.ensembl.org"

    def __init__(self):
        pass

    def _make_request(self, endpoint: str, content_type: str = "application/json") -> Optional[Any]:
        """Make a request to the Ensembl REST API."""
        url = f"{self.BASE_URL}{endpoint}"

        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Content-Type": content_type,
                    "Accept": content_type
                }
            )
            with urllib.request.urlopen(req, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as e:
            return {"error": str(e)}

    def lookup_gene(self, gene_symbol: str, species: str = "human") -> Optional[Dict[str, Any]]:
        """Look up a gene by symbol.

        Args:
            gene_symbol: Gene symbol (e.g., "BRCA1", "TP53")
            species: Species name (default: "human")

        Returns:
            Gene information dictionary or None
        """
        endpoint = f"/lookup/symbol/{species}/{gene_symbol}?expand=1"
        result = self._make_request(endpoint)

        if result and "error" not in result:
            return result
        return None

    def lookup_by_id(self, ensembl_id: str) -> Optional[Dict[str, Any]]:
        """Look up an Ensembl ID.

        Args:
            ensembl_id: Ensembl ID (e.g., "ENSG00000012048")

        Returns:
            Gene/transcript information dictionary or None
        """
        endpoint = f"/lookup/id/{ensembl_id}?expand=1"
        result = self._make_request(endpoint)

        if result and "error" not in result:
            return result
        return None

    def get_sequence(self, ensembl_id: str, seq_type: str = "genomic") -> Optional[Dict[str, Any]]:
        """Get sequence for an Ensembl ID.

        Args:
            ensembl_id: Ensembl ID
            seq_type: Sequence type (genomic, cds, cdna, protein)

        Returns:
            Sequence dictionary or None
        """
        endpoint = f"/sequence/id/{ensembl_id}?type={seq_type}"
        result = self._make_request(endpoint)

        if result and "error" not in result:
            return result
        return None

    def get_region(self, species: str, region: str) -> Optional[List[Dict[str, Any]]]:
        """Get features in a genomic region.

        Args:
            species: Species name (e.g., "human")
            region: Genomic region (e.g., "17:43044295-43125483")

        Returns:
            List of features or None
        """
        endpoint = f"/overlap/region/{species}/{region}?feature=gene"
        result = self._make_request(endpoint)

        if result and isinstance(result, list):
            return result
        return None

    def vep_variant(self, species: str, variant: str) -> Optional[List[Dict[str, Any]]]:
        """Run VEP (Variant Effect Predictor) on a variant.

        Args:
            species: Species name (e.g., "human")
            variant: Variant in format "chr:pos:ref:alt" or HGVS notation

        Returns:
            VEP results or None
        """
        # Parse variant format
        if "/" in variant:
            # HGVS format
            endpoint = f"/vep/{species}/hgvs/{urllib.parse.quote(variant)}"
        else:
            # chr:pos:ref:alt format
            parts = variant.replace(":", " ").replace("-", " ").split()
            if len(parts) >= 4:
                chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
                endpoint = f"/vep/{species}/region/{chrom}:{pos}:{pos}/{alt}"
            else:
                return None

        result = self._make_request(endpoint)

        if result and isinstance(result, list):
            return result
        return None

    def search_genes(self, query: str, species: str = "human") -> Optional[List[Dict[str, Any]]]:
        """Search for genes by name or description.

        Args:
            query: Search query
            species: Species name

        Returns:
            List of matching genes or None
        """
        # Use xrefs endpoint for search
        endpoint = f"/xrefs/symbol/{species}/{query}"
        result = self._make_request(endpoint)

        if result and isinstance(result, list):
            # Get full details for each result
            genes = []
            for xref in result[:5]:  # Limit to 5 results
                gene_id = xref.get("id")
                if gene_id:
                    gene_info = self.lookup_by_id(gene_id)
                    if gene_info:
                        genes.append(gene_info)
            return genes
        return None


def format_gene_info(gene: Dict[str, Any]) -> str:
    """Format gene information for display."""
    if not gene:
        return "Gene not found."

    # Basic info
    output = f"""
## Gene: {gene.get('display_name', gene.get('id', 'Unknown'))}

- **Ensembl ID**: {gene.get('id', 'N/A')}
- **Biotype**: {gene.get('biotype', 'N/A')}
- **Description**: {gene.get('description', 'N/A')}
- **Species**: {gene.get('species', 'N/A')}

### Location
- **Chromosome**: {gene.get('seq_region_name', 'N/A')}
- **Start**: {gene.get('start', 'N/A')}
- **End**: {gene.get('end', 'N/A')}
- **Strand**: {'+' if gene.get('strand', 1) > 0 else '-'}
- **Assembly**: {gene.get('assembly_name', 'N/A')}
"""

    # Transcripts
    transcripts = gene.get('Transcript', [])
    if transcripts:
        output += f"\n### Transcripts ({len(transcripts)} total)\n"
        for t in transcripts[:5]:  # Show first 5
            output += f"- **{t.get('display_name', t.get('id'))}** ({t.get('biotype', 'N/A')})\n"

    return output


def format_vep_results(results: List[Dict[str, Any]]) -> str:
    """Format VEP results for display."""
    if not results:
        return "No VEP results found."

    output = []
    for result in results:
        variant_output = f"""
## Variant: {result.get('input', 'Unknown')}

- **Most Severe Consequence**: {result.get('most_severe_consequence', 'N/A')}
- **Colocated Variants**: {len(result.get('colocated_variants', []))}
"""

        # Transcript consequences
        consequences = result.get('transcript_consequences', [])
        if consequences:
            variant_output += "\n### Transcript Consequences\n"
            for c in consequences[:5]:  # Show first 5
                variant_output += f"""
- **Gene**: {c.get('gene_symbol', 'N/A')} ({c.get('gene_id', 'N/A')})
  - Consequence: {', '.join(c.get('consequence_terms', ['N/A']))}
  - Impact: {c.get('impact', 'N/A')}
  - Amino acids: {c.get('amino_acids', 'N/A')}
  - Codons: {c.get('codons', 'N/A')}
"""

        output.append(variant_output)

    return "\n---\n".join(output)


def format_region_results(features: List[Dict[str, Any]]) -> str:
    """Format region query results for display."""
    if not features:
        return "No features found in region."

    output = [f"## Found {len(features)} genes in region\n"]

    for f in features[:10]:  # Show first 10
        output.append(f"""
- **{f.get('external_name', f.get('id', 'Unknown'))}**
  - ID: {f.get('id', 'N/A')}
  - Biotype: {f.get('biotype', 'N/A')}
  - Location: {f.get('start', 'N/A')}-{f.get('end', 'N/A')}
""")

    return "\n".join(output)


def main():
    """Main function for MCP server."""
    client = EnsemblClient()

    # Read from stdin for MCP protocol
    for line in sys.stdin:
        try:
            request = json.loads(line.strip())

            method = request.get("method", "")
            params = request.get("params", {})

            if method == "lookup_gene":
                gene_symbol = params.get("gene", params.get("symbol", ""))
                species = params.get("species", "human")

                gene = client.lookup_gene(gene_symbol, species)
                result = format_gene_info(gene)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "found": gene is not None
                    }
                }

            elif method == "lookup_id":
                ensembl_id = params.get("id", "")

                gene = client.lookup_by_id(ensembl_id)
                result = format_gene_info(gene)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "found": gene is not None
                    }
                }

            elif method == "get_region":
                species = params.get("species", "human")
                region = params.get("region", "")

                features = client.get_region(species, region)
                result = format_region_results(features) if features else "No features found."

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(features) if features else 0
                    }
                }

            elif method == "vep":
                species = params.get("species", "human")
                variant = params.get("variant", "")

                vep_results = client.vep_variant(species, variant)
                result = format_vep_results(vep_results) if vep_results else "VEP analysis failed."

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(vep_results) if vep_results else 0
                    }
                }

            elif method == "search":
                query = params.get("query", "")
                species = params.get("species", "human")

                genes = client.search_genes(query, species)
                if genes:
                    result = "\n---\n".join([format_gene_info(g) for g in genes])
                else:
                    result = "No genes found."

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(genes) if genes else 0
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
        client = EnsemblClient()
        query = " ".join(sys.argv[2:]) if len(sys.argv) > 2 else "BRCA1"
        print(f"Looking up gene: {query}")
        gene = client.lookup_gene(query)
        if gene:
            print(format_gene_info(gene))
        else:
            print(f"Gene '{query}' not found.")
    else:
        main()
