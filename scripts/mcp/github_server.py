#!/usr/bin/env python3
"""
GitHub MCP Server for InterSubMod Project

This server provides GitHub search capabilities via the Model Context Protocol (MCP).
It uses the GitHub REST API to search repositories, code, and issues.

Usage:
    The server will automatically load GitHub token from:
    1. Environment variable GITHUB_TOKEN
    2. .env file in the project root

    python github_server.py

    Test mode:
    python github_server.py --test "methylation analysis"

Environment Variables:
    GITHUB_TOKEN: Your GitHub personal access token for higher rate limits
"""

import json
import os
import sys
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Any, Dict, List, Optional


def load_env_file():
    """Load environment variables from .env file."""
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


class GitHubClient:
    """Client for interacting with GitHub REST API."""

    BASE_URL = "https://api.github.com"

    def __init__(self, token: Optional[str] = None):
        self.token = token or os.environ.get("GITHUB_TOKEN", "")

    def _make_request(self, endpoint: str) -> Optional[Dict]:
        """Make a request to the GitHub API."""
        url = f"{self.BASE_URL}{endpoint}"

        headers = {
            "Accept": "application/vnd.github.v3+json",
            "User-Agent": "InterSubMod-MCP/1.0"
        }

        if self.token:
            headers["Authorization"] = f"token {self.token}"

        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            return {"error": f"HTTP {e.code}: {e.reason}"}
        except Exception as e:
            return {"error": str(e)}

    def search_repositories(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """Search GitHub repositories.

        Args:
            query: Search query
            max_results: Maximum number of results

        Returns:
            List of repository dictionaries
        """
        encoded_query = urllib.parse.quote(query)
        endpoint = f"/search/repositories?q={encoded_query}&sort=stars&per_page={max_results}"

        result = self._make_request(endpoint)

        if result and "items" in result:
            return result["items"]
        return []

    def search_code(self, query: str, language: Optional[str] = None, max_results: int = 10) -> List[Dict[str, Any]]:
        """Search code on GitHub.

        Args:
            query: Search query
            language: Programming language filter
            max_results: Maximum number of results

        Returns:
            List of code result dictionaries
        """
        search_query = query
        if language:
            search_query += f"+language:{language}"

        encoded_query = urllib.parse.quote(search_query)
        endpoint = f"/search/code?q={encoded_query}&per_page={max_results}"

        result = self._make_request(endpoint)

        if result and "items" in result:
            return result["items"]
        return []

    def search_issues(self, query: str, repo: Optional[str] = None, max_results: int = 10) -> List[Dict[str, Any]]:
        """Search GitHub issues and pull requests.

        Args:
            query: Search query
            repo: Repository filter (e.g., "owner/repo")
            max_results: Maximum number of results

        Returns:
            List of issue/PR dictionaries
        """
        search_query = query
        if repo:
            search_query += f"+repo:{repo}"

        encoded_query = urllib.parse.quote(search_query)
        endpoint = f"/search/issues?q={encoded_query}&per_page={max_results}"

        result = self._make_request(endpoint)

        if result and "items" in result:
            return result["items"]
        return []

    def get_repo_readme(self, owner: str, repo: str) -> Optional[str]:
        """Get README content for a repository.

        Args:
            owner: Repository owner
            repo: Repository name

        Returns:
            README content or None
        """
        endpoint = f"/repos/{owner}/{repo}/readme"

        result = self._make_request(endpoint)

        if result and "content" in result:
            import base64
            try:
                return base64.b64decode(result["content"]).decode("utf-8")
            except Exception:
                return None
        return None

    def get_repo_info(self, owner: str, repo: str) -> Optional[Dict[str, Any]]:
        """Get repository information.

        Args:
            owner: Repository owner
            repo: Repository name

        Returns:
            Repository information dictionary or None
        """
        endpoint = f"/repos/{owner}/{repo}"

        result = self._make_request(endpoint)

        if result and "error" not in result:
            return result
        return None


def format_repo_results(repos: List[Dict[str, Any]]) -> str:
    """Format repository search results for display."""
    if not repos:
        return "No repositories found."

    output = []
    for i, repo in enumerate(repos, 1):
        description = repo.get("description", "No description")
        if description and len(description) > 200:
            description = description[:200] + "..."

        output.append(f"""
## {i}. {repo.get('full_name', 'Unknown')}

- **Stars**: {repo.get('stargazers_count', 0):,}
- **Forks**: {repo.get('forks_count', 0):,}
- **Language**: {repo.get('language', 'N/A')}
- **Last Updated**: {repo.get('updated_at', 'N/A')[:10] if repo.get('updated_at') else 'N/A'}
- **URL**: {repo.get('html_url', 'N/A')}

**Description**: {description}
""")

    return "\n---\n".join(output)


def format_code_results(results: List[Dict[str, Any]]) -> str:
    """Format code search results for display."""
    if not results:
        return "No code results found."

    output = []
    for i, item in enumerate(results, 1):
        repo = item.get("repository", {})

        output.append(f"""
## {i}. {item.get('name', 'Unknown')}

- **Repository**: {repo.get('full_name', 'N/A')}
- **Path**: {item.get('path', 'N/A')}
- **URL**: {item.get('html_url', 'N/A')}
""")

    return "\n---\n".join(output)


def format_issue_results(issues: List[Dict[str, Any]]) -> str:
    """Format issue/PR search results for display."""
    if not issues:
        return "No issues found."

    output = []
    for i, issue in enumerate(issues, 1):
        title = issue.get("title", "No title")
        if len(title) > 100:
            title = title[:100] + "..."

        body = issue.get("body", "") or ""
        if len(body) > 200:
            body = body[:200] + "..."

        is_pr = "pull_request" in issue
        issue_type = "PR" if is_pr else "Issue"

        output.append(f"""
## {i}. [{issue_type}] {title}

- **Repository**: {issue.get('repository_url', 'N/A').split('repos/')[-1] if issue.get('repository_url') else 'N/A'}
- **State**: {issue.get('state', 'N/A')}
- **Created**: {issue.get('created_at', 'N/A')[:10] if issue.get('created_at') else 'N/A'}
- **URL**: {issue.get('html_url', 'N/A')}

{body if body else 'No description'}
""")

    return "\n---\n".join(output)


def format_repo_info(repo: Dict[str, Any]) -> str:
    """Format repository information for display."""
    if not repo:
        return "Repository not found."

    description = repo.get("description", "No description")

    return f"""
## Repository: {repo.get('full_name', 'Unknown')}

- **Stars**: {repo.get('stargazers_count', 0):,}
- **Forks**: {repo.get('forks_count', 0):,}
- **Watchers**: {repo.get('watchers_count', 0):,}
- **Open Issues**: {repo.get('open_issues_count', 0):,}
- **Language**: {repo.get('language', 'N/A')}
- **License**: {repo.get('license', {}).get('name', 'N/A') if repo.get('license') else 'N/A'}
- **Created**: {repo.get('created_at', 'N/A')[:10] if repo.get('created_at') else 'N/A'}
- **Last Updated**: {repo.get('updated_at', 'N/A')[:10] if repo.get('updated_at') else 'N/A'}
- **URL**: {repo.get('html_url', 'N/A')}

**Description**: {description}

**Topics**: {', '.join(repo.get('topics', [])) or 'None'}
"""


def main():
    """Main function for MCP server."""
    client = GitHubClient()

    # Read from stdin for MCP protocol
    for line in sys.stdin:
        try:
            request = json.loads(line.strip())

            method = request.get("method", "")
            params = request.get("params", {})

            if method == "search_repos":
                query = params.get("query", "")
                max_results = params.get("max_results", 10)

                repos = client.search_repositories(query, max_results)
                result = format_repo_results(repos)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(repos)
                    }
                }

            elif method == "search_code":
                query = params.get("query", "")
                language = params.get("language")
                max_results = params.get("max_results", 10)

                results = client.search_code(query, language, max_results)
                result = format_code_results(results)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(results)
                    }
                }

            elif method == "search_issues":
                query = params.get("query", "")
                repo = params.get("repo")
                max_results = params.get("max_results", 10)

                issues = client.search_issues(query, repo, max_results)
                result = format_issue_results(issues)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "count": len(issues)
                    }
                }

            elif method == "get_repo":
                owner = params.get("owner", "")
                repo_name = params.get("repo", "")

                # Handle "owner/repo" format
                if "/" in owner and not repo_name:
                    owner, repo_name = owner.split("/", 1)

                repo = client.get_repo_info(owner, repo_name)
                result = format_repo_info(repo)

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "found": repo is not None
                    }
                }

            elif method == "get_readme":
                owner = params.get("owner", "")
                repo_name = params.get("repo", "")

                # Handle "owner/repo" format
                if "/" in owner and not repo_name:
                    owner, repo_name = owner.split("/", 1)

                readme = client.get_repo_readme(owner, repo_name)
                result = readme if readme else "README not found."

                response = {
                    "jsonrpc": "2.0",
                    "id": request.get("id"),
                    "result": {
                        "content": result,
                        "found": readme is not None
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
        client = GitHubClient()
        query = " ".join(sys.argv[2:]) if len(sys.argv) > 2 else "methylation analysis"
        print(f"Searching GitHub for: {query}")
        repos = client.search_repositories(query, 5)
        print(f"Found {len(repos)} repositories")
        print(format_repo_results(repos))
    else:
        main()
