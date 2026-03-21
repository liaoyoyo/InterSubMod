#!/usr/bin/env python3
"""
Shared Planning Bridge Script

A cross-platform utility for managing task planning files that can be used by:
- Claude Code (via skill system)
- Google Antigravity (via custom agent or direct invocation)

This script provides a unified interface for creating and updating planning files,
ensuring consistency across different AI development platforms.

Usage:
    python bridge.py init <task_name> [--type TYPE]
    python bridge.py update <task_name> --progress "message"
    python bridge.py update <task_name> --finding "title" "content"
    python bridge.py update <task_name> --complete-task "task_id"
    python bridge.py list
    python bridge.py stats [task_name]
    python bridge.py validate <task_name>
    python bridge.py unlock <task_name>
"""

import os
import sys
import json
import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re


class SharedPlanningBridge:
    """Bridge for cross-platform task planning."""

    def __init__(self, project_root: Optional[str] = None):
        """Initialize the bridge with project root directory."""
        if project_root is None:
            # Auto-detect project root (look for .git directory)
            current = Path.cwd()
            while current != current.parent:
                if (current / '.git').exists():
                    project_root = str(current)
                    break
                current = current.parent
            else:
                project_root = str(Path.cwd())

        self.project_root = Path(project_root)
        self.config_path = self.project_root / '.claude/skills/shared-planning/config.json'
        self.config = self._load_config()
        self.platform = self._detect_platform()

    def _load_config(self) -> Dict:
        """Load configuration from config.json."""
        if not self.config_path.exists():
            # Use default config
            return {
                "paths": {
                    "base_dir": "docs/shared_planning",
                    "task_plan_file": "task_plan.md",
                    "findings_file": "findings.md",
                    "progress_file": "progress.md"
                },
                "platforms": {
                    "claude_code": {"identifier": "Claude Code"},
                    "antigravity": {"identifier": "Antigravity"}
                }
            }

        with open(self.config_path, 'r', encoding='utf-8') as f:
            return json.load(f)

    def _detect_platform(self) -> str:
        """Detect which platform is running this script."""
        platforms = self.config.get('platforms', {})

        # Check environment variables
        if os.environ.get('CLAUDE_CODE'):
            return platforms.get('claude_code', {}).get('identifier', 'Claude Code')
        elif os.environ.get('ANTIGRAVITY_SESSION'):
            return platforms.get('antigravity', {}).get('identifier', 'Antigravity')

        # If no env var, check for interactive session markers
        # This is a heuristic and may need adjustment
        if 'ANTHROPIC' in os.environ or 'CLAUDE' in os.environ:
            return 'Claude Code'

        return 'Unknown'

    def _get_task_dir(self, task_name: str) -> Path:
        """Get the directory path for a task."""
        base_dir = self.config['paths']['base_dir']
        return self.project_root / base_dir / self._sanitize_task_name(task_name)

    def _sanitize_task_name(self, task_name: str) -> str:
        """Sanitize task name for use as directory name."""
        # Replace spaces with underscores, remove special chars
        sanitized = re.sub(r'[^\w\s-]', '', task_name.lower())
        sanitized = re.sub(r'[-\s]+', '_', sanitized)
        return sanitized.strip('_')

    def _get_timestamp(self) -> str:
        """Get current timestamp in standard format."""
        return datetime.now().strftime('%Y-%m-%d %H:%M')

    def _get_date(self) -> str:
        """Get current date."""
        return datetime.now().strftime('%Y-%m-%d')

    def _get_day_of_week(self) -> str:
        """Get current day of week."""
        days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
        return days[datetime.now().weekday()]

    def init_task(self, task_name: str, task_type: str = 'feature',
                  objective: str = '') -> Tuple[bool, str]:
        """
        Initialize a new task with planning files.

        Args:
            task_name: Name of the task
            task_type: Type of task (research, feature, bugfix, etc.)
            objective: Optional objective description

        Returns:
            Tuple of (success: bool, message: str)
        """
        task_dir = self._get_task_dir(task_name)

        # Check if task already exists
        if task_dir.exists():
            return False, f"Task '{task_name}' already exists at {task_dir}"

        # Create task directory
        task_dir.mkdir(parents=True, exist_ok=True)

        # Create planning files
        timestamp = self._get_timestamp()
        date = self._get_date()

        # Create task_plan.md
        task_plan_content = self._generate_task_plan_template(
            task_name, task_type, timestamp, objective
        )
        task_plan_path = task_dir / self.config['paths']['task_plan_file']
        task_plan_path.write_text(task_plan_content, encoding='utf-8')

        # Create findings.md
        findings_content = self._generate_findings_template(task_name, date)
        findings_path = task_dir / self.config['paths']['findings_file']
        findings_path.write_text(findings_content, encoding='utf-8')

        # Create progress.md
        progress_content = self._generate_progress_template(task_name, date)
        progress_path = task_dir / self.config['paths']['progress_file']
        progress_path.write_text(progress_content, encoding='utf-8')

        return True, f"Successfully initialized task '{task_name}' at {task_dir}"

    def _generate_task_plan_template(self, task_name: str, task_type: str,
                                     timestamp: str, objective: str = '') -> str:
        """Generate task_plan.md template."""
        obj_text = objective if objective else "[Clear, concise description of the task goal]"

        return f"""# Task Plan: {task_name}

**Created**: {timestamp} (via {self.platform})
**Updated**: {timestamp} (via {self.platform})
**Status**: In Progress
**Type**: {task_type}

## Objective

{obj_text}

## Background

[Context and prerequisites for this task]

## Tasks

### Phase 1: Planning & Research
- [ ] Define requirements
  - Assigned to: {self.platform}
  - Priority: High
- [ ] Research existing solutions
  - Assigned to: {self.platform}
  - Priority: Medium

### Phase 2: Implementation
- [ ] Implement core functionality
  - Assigned to: [Platform]
  - Priority: High

### Phase 3: Testing & Validation
- [ ] Write tests
  - Assigned to: [Platform]
  - Priority: High
- [ ] Validate results
  - Assigned to: [Platform]
  - Priority: High

## Dependencies

- **Depends on**: [List external dependencies or prerequisite tasks]
- **Blocks**: [List tasks that depend on this]

## Resources

- **Documentation**: [Links to relevant documentation]
- **Related Files**:
  - `path/to/file` - Description
- **References**: [Links to papers, articles, etc.]

## Notes

[Any additional notes or considerations]
"""

    def _generate_findings_template(self, task_name: str, date: str) -> str:
        """Generate findings.md template."""
        return f"""# Findings: {task_name}

**Created**: {date} (via {self.platform})
**Last Updated**: {date} (via {self.platform})

## Key Discoveries

_Document important findings here as they are discovered._

### Example Finding (Remove this section)

**Context**: What question were we trying to answer?

**Discovery**: What did we find?

**Evidence**:
- Code location: `src/file.cpp:123`
- Test results: [Summary or link]
- Data analysis: [Link to charts/tables]

**Implications**:
- How does this affect our approach?
- What decisions can we make based on this?

**Action Items**:
- [ ] Follow-up action 1
- [ ] Follow-up action 2

---

## Research Notes

### Topic 1: [Topic Name]

- Note 1
- Note 2

### References

[1] Citation format: Author (Year). "Title". Source. URL
"""

    def _generate_progress_template(self, task_name: str, date: str) -> str:
        """Generate progress.md template."""
        day = self._get_day_of_week()
        return f"""# Progress Log: {task_name}

**Started**: {date}
**Current Status**: 0% complete | In Progress
**Last Updated**: {date} (via {self.platform})

---

## {date} {day}

**Session Info**:
- Platform: {self.platform}
- Focus area: Task initialization

### ✅ Completed Today
- **Task initialization**: Created planning files
  - Files created: task_plan.md, findings.md, progress.md

### 🔄 In Progress
- **Planning**: Defining detailed task breakdown

### 📋 Planned Next
- **Research**: [Next steps]

### 💡 Key Insights
- Task initialized and ready for work

### ⏱️ Time Allocation
- Planning: 0.5 hours

---
"""

    def update_progress(self, task_name: str, message: str) -> Tuple[bool, str]:
        """Add a progress entry to progress.md."""
        task_dir = self._get_task_dir(task_name)
        progress_path = task_dir / self.config['paths']['progress_file']

        if not progress_path.exists():
            return False, f"Task '{task_name}' not found. Initialize it first with 'init' command."

        # Read existing content
        content = progress_path.read_text(encoding='utf-8')

        # Update the "Last Updated" line
        timestamp = self._get_timestamp()
        content = re.sub(
            r'\*\*Last Updated\*\*: .*',
            f'**Last Updated**: {timestamp} (via {self.platform})',
            content
        )

        # Add new progress entry at the top (after the header)
        date = self._get_date()
        day = self._get_day_of_week()

        new_entry = f"""## {date} {day}

**Session Info**:
- Platform: {self.platform}
- Update: {message}

### ✅ Completed Today
- {message}

---

"""

        # Insert after the first "---" separator (after header)
        parts = content.split('---', 2)
        if len(parts) >= 2:
            content = parts[0] + '---' + parts[1] + '---\n\n' + new_entry + parts[2]
        else:
            content += '\n\n' + new_entry

        # Write back
        progress_path.write_text(content, encoding='utf-8')

        return True, f"Progress updated for task '{task_name}'"

    def add_finding(self, task_name: str, title: str, context: str = '',
                    discovery: str = '', evidence: List[str] = None,
                    implications: List[str] = None) -> Tuple[bool, str]:
        """Add a finding to findings.md."""
        task_dir = self._get_task_dir(task_name)
        findings_path = task_dir / self.config['paths']['findings_file']

        if not findings_path.exists():
            return False, f"Task '{task_name}' not found."

        # Read existing content
        content = findings_path.read_text(encoding='utf-8')

        # Update timestamp
        date = self._get_date()
        content = re.sub(
            r'\*\*Last Updated\*\*: .*',
            f'**Last Updated**: {date} (via {self.platform})',
            content
        )

        # Create new finding entry
        evidence_list = '\n'.join([f'- {e}' for e in (evidence or [])])
        implications_list = '\n'.join([f'- {i}' for i in (implications or [])])

        new_finding = f"""### {date}: {title}

**Context**: {context or '[Describe the context]'}

**Discovery**: {discovery or '[Describe what was found]'}

**Evidence**:
{evidence_list or '- [Add evidence]'}

**Implications**:
{implications_list or '- [Describe implications]'}

**Action Items**:
- [ ] [Add action items]

---

"""

        # Insert after "## Key Discoveries" section
        if '## Key Discoveries' in content:
            parts = content.split('## Key Discoveries', 1)
            # Find the end of the section header (next line)
            after_header = parts[1].split('\n', 2)
            if len(after_header) >= 2:
                content = (parts[0] + '## Key Discoveries' + '\n' +
                          after_header[1] + '\n\n' + new_finding +
                          after_header[2] if len(after_header) > 2 else '')
            else:
                content = parts[0] + '## Key Discoveries\n\n' + new_finding + parts[1]
        else:
            # Append to end
            content += '\n\n' + new_finding

        # Write back
        findings_path.write_text(content, encoding='utf-8')

        return True, f"Finding '{title}' added to task '{task_name}'"

    def complete_task_item(self, task_name: str, task_pattern: str) -> Tuple[bool, str]:
        """Mark a task item as completed in task_plan.md."""
        task_dir = self._get_task_dir(task_name)
        task_plan_path = task_dir / self.config['paths']['task_plan_file']

        if not task_plan_path.exists():
            return False, f"Task '{task_name}' not found."

        # Read content
        content = task_plan_path.read_text(encoding='utf-8')

        # Find and update the task item
        # Look for "- [ ] {task_pattern}" and replace with "- [x] {task_pattern}"
        pattern = re.compile(
            rf'^(\s*)- \[ \] (.*)' + re.escape(task_pattern) + r'(.*)',
            re.MULTILINE
        )

        if not pattern.search(content):
            return False, f"Task item matching '{task_pattern}' not found"

        # Replace with completed checkbox
        date = self._get_date()
        content = pattern.sub(
            rf'\1- [x] \2{re.escape(task_pattern)}\3\n\1  - Completed: {date} (via {self.platform})',
            content
        )

        # Update timestamp
        timestamp = self._get_timestamp()
        content = re.sub(
            r'\*\*Updated\*\*: .*',
            f'**Updated**: {timestamp} (via {self.platform})',
            content
        )

        # Write back
        task_plan_path.write_text(content, encoding='utf-8')

        return True, f"Task item '{task_pattern}' marked as completed"

    def list_tasks(self) -> List[Dict]:
        """List all tasks."""
        base_dir = self.project_root / self.config['paths']['base_dir']

        if not base_dir.exists():
            return []

        tasks = []
        for task_dir in base_dir.iterdir():
            if task_dir.is_dir():
                task_plan_path = task_dir / self.config['paths']['task_plan_file']
                if task_plan_path.exists():
                    # Parse task info
                    content = task_plan_path.read_text(encoding='utf-8')

                    # Extract status
                    status_match = re.search(r'\*\*Status\*\*: (.+)', content)
                    status = status_match.group(1) if status_match else 'Unknown'

                    # Extract updated time
                    updated_match = re.search(r'\*\*Updated\*\*: (.+)', content)
                    updated = updated_match.group(1) if updated_match else 'Unknown'

                    tasks.append({
                        'name': task_dir.name,
                        'path': str(task_dir),
                        'status': status,
                        'updated': updated
                    })

        return tasks

    def get_task_stats(self, task_name: Optional[str] = None) -> Dict:
        """Get statistics for a task or all tasks."""
        if task_name:
            # Stats for specific task
            task_dir = self._get_task_dir(task_name)
            task_plan_path = task_dir / self.config['paths']['task_plan_file']

            if not task_plan_path.exists():
                return {'error': f"Task '{task_name}' not found"}

            content = task_plan_path.read_text(encoding='utf-8')

            # Count tasks
            total_tasks = len(re.findall(r'^- \[[ x]\]', content, re.MULTILINE))
            completed_tasks = len(re.findall(r'^- \[x\]', content, re.MULTILINE))

            return {
                'task_name': task_name,
                'total_tasks': total_tasks,
                'completed_tasks': completed_tasks,
                'completion_rate': f"{completed_tasks/total_tasks*100:.1f}%" if total_tasks > 0 else "0%"
            }
        else:
            # Stats for all tasks
            tasks = self.list_tasks()
            total_tasks = len(tasks)

            status_counts = {}
            for task in tasks:
                status = task['status']
                status_counts[status] = status_counts.get(status, 0) + 1

            return {
                'total_tasks': total_tasks,
                'by_status': status_counts,
                'tasks': tasks
            }

    def validate_task(self, task_name: str) -> Tuple[bool, List[str]]:
        """Validate task files for consistency."""
        task_dir = self._get_task_dir(task_name)
        issues = []

        # Check if task exists
        if not task_dir.exists():
            return False, [f"Task '{task_name}' does not exist"]

        # Check required files
        required_files = [
            self.config['paths']['task_plan_file'],
            self.config['paths']['findings_file'],
            self.config['paths']['progress_file']
        ]

        for file in required_files:
            file_path = task_dir / file
            if not file_path.exists():
                issues.append(f"Missing required file: {file}")
            else:
                # Validate file format
                content = file_path.read_text(encoding='utf-8')

                # Check for required headers
                if file == 'task_plan.md':
                    required = ['Created', 'Updated', 'Status', 'Objective']
                elif file == 'findings.md':
                    required = ['Created', 'Last Updated']
                else:  # progress.md
                    required = ['Started', 'Current Status']

                for req in required:
                    if f'**{req}**' not in content:
                        issues.append(f"{file}: Missing required field '{req}'")

        return len(issues) == 0, issues


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Shared Planning Bridge - Cross-platform task planning utility'
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Init command
    init_parser = subparsers.add_parser('init', help='Initialize a new task')
    init_parser.add_argument('task_name', help='Name of the task')
    init_parser.add_argument('--type', default='feature',
                            choices=['research', 'feature', 'bugfix', 'refactor', 'experiment', 'documentation'],
                            help='Type of task')
    init_parser.add_argument('--objective', default='', help='Task objective')

    # Update command
    update_parser = subparsers.add_parser('update', help='Update task files')
    update_parser.add_argument('task_name', help='Name of the task')
    update_parser.add_argument('--progress', help='Add progress message')
    update_parser.add_argument('--finding-title', help='Finding title')
    update_parser.add_argument('--finding-context', default='', help='Finding context')
    update_parser.add_argument('--finding-discovery', default='', help='Finding discovery')
    update_parser.add_argument('--complete-task', help='Mark task as complete (pattern to match)')

    # List command
    list_parser = subparsers.add_parser('list', help='List all tasks')

    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Show task statistics')
    stats_parser.add_argument('task_name', nargs='?', help='Optional: specific task name')

    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate task files')
    validate_parser.add_argument('task_name', help='Name of the task')

    args = parser.parse_args()

    # Initialize bridge
    bridge = SharedPlanningBridge()

    # Execute command
    if args.command == 'init':
        success, message = bridge.init_task(args.task_name, args.type, args.objective)
        print(message)
        sys.exit(0 if success else 1)

    elif args.command == 'update':
        if args.progress:
            success, message = bridge.update_progress(args.task_name, args.progress)
            print(message)
            sys.exit(0 if success else 1)
        elif args.finding_title:
            success, message = bridge.add_finding(
                args.task_name, args.finding_title,
                context=args.finding_context,
                discovery=args.finding_discovery
            )
            print(message)
            sys.exit(0 if success else 1)
        elif args.complete_task:
            success, message = bridge.complete_task_item(args.task_name, args.complete_task)
            print(message)
            sys.exit(0 if success else 1)
        else:
            print("Error: --progress, --finding-title, or --complete-task required for update command")
            sys.exit(1)

    elif args.command == 'list':
        tasks = bridge.list_tasks()
        if not tasks:
            print("No tasks found")
        else:
            print(f"\nFound {len(tasks)} task(s):\n")
            for task in tasks:
                print(f"  • {task['name']}")
                print(f"    Status: {task['status']}")
                print(f"    Updated: {task['updated']}")
                print(f"    Path: {task['path']}\n")

    elif args.command == 'stats':
        stats = bridge.get_task_stats(args.task_name)
        if 'error' in stats:
            print(f"Error: {stats['error']}")
            sys.exit(1)
        else:
            print("\nTask Statistics:\n")
            print(json.dumps(stats, indent=2))

    elif args.command == 'validate':
        valid, issues = bridge.validate_task(args.task_name)
        if valid:
            print(f"✓ Task '{args.task_name}' is valid")
            sys.exit(0)
        else:
            print(f"✗ Task '{args.task_name}' has issues:\n")
            for issue in issues:
                print(f"  • {issue}")
            sys.exit(1)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
