#!/usr/bin/env python3
import argparse
import glob
import os
import shutil
import subprocess
import sys

# --- Configuration ---
GEMINI_CMD = "gemini"  # The executable name for the Gemini CLI

# The review prompt containing To-Dos for a high-quality review
REVIEW_PROMPT = (
    "You are an expert code reviewer. Please review the provided code file according to the following To-Dos: "
    "1. **Correctness**: Identify any logical errors, bugs, or race conditions. "
    "2. **Style & Conventions**: Ensure the code follows project conventions (PEP8 for Python) and idiomatic patterns. "
    "3. **Performance**: Point out any inefficient algorithms or unnecessary resource usage. "
    "4. **Readability**: Suggest improvements for variable naming, function size, and clarity. "
    "5. **Safety**: Check for security vulnerabilities or improper error handling. "
    "6. **Documentation**: Verify that complex logic is adequately explained (focusing on 'why', not 'what'). "
    "Write all findings to the folder output/*filename*_review_findings.md as todos"
)

def get_files(args):
    """Collects a list of files based on arguments (recursive or glob patterns)."""
    files = []
    if args.recursive:
        # Recursive search for all files, ignoring common non-source directories
        for root, _, filenames in os.walk("."):
            if any(x in root for x in [".git", "__pycache__", ".venv", ".vscode", "output"]):
                continue
            for f in filenames:
                files.append(os.path.join(root, f))
    elif args.files:
        for pattern in args.files:
            # Shells usually expand globs, but we handle explicit glob strings just in case
            matched = glob.glob(pattern)
            if matched:
                files.extend(matched)
            else:
                # If no glob match, assume it's a specific file path
                files.append(pattern)

    # Remove duplicates and sort
    return sorted(list(set(files)))

def review_file(file_path):
    """Starts a Gemini CLI instance for the specified file."""
    if not os.path.isfile(file_path):
        print(f"Skipping '{file_path}': Not a file.")
        return

    print(f"\\n{'='*60}")
    print(f"REVIEWING: {file_path}")
    print(f"{ '='*60}")

    # We append the specific file path to the prompt so the agent knows what to look at.
    # The agent can then use its tools (read_file) to inspect the content.
    full_prompt = f"{REVIEW_PROMPT} Please review the file: {file_path} IMPORTANT: Provide the review as text output only. Do NOT attempt to create or modify any files."

    try:
        # Resolve the executable path (handles .cmd/.bat on Windows)
        executable = shutil.which(GEMINI_CMD) or GEMINI_CMD

        print(f"Launching {executable} for {file_path}...")
        # Run in one-shot mode, capturing output.
        # Explicitly set encoding to utf-8 to handle emoji/special chars on Windows
        result = subprocess.run(
            [executable, full_prompt],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace"
        )

        if result.returncode == 0:
            # Construct path inside the output/ directory
            # We use normpath and lstrip to handle potential absolute paths safely
            rel_path = os.path.normpath(file_path).lstrip(os.sep)
            output_filename = os.path.join("output", f"{rel_path}.review.md")

            # Ensure the target subdirectory exists
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)

            if result.stdout:
                with open(output_filename, "w", encoding="utf-8") as f:
                    f.write(result.stdout)
                print(f"Review saved to: {output_filename}")
            else:
                print("Warning: Review successful but no output captured.")
        else:
            print(f"Error reviewing {file_path}:")
            print(result.stderr)

    except FileNotFoundError:
        print(f"Error: Command '{GEMINI_CMD}' not found. Please ensure it is installed and in your PATH.")
    except KeyboardInterrupt:
        print("\nSession interrupted.")

def main():
    parser = argparse.ArgumentParser(description="Iteratively review files using the Gemini CLI.")        
    parser.add_argument("-r", "--recursive", action="store_true", help="Recursively review all files in the current directory.")
    parser.add_argument("files", nargs="*", help="Files or glob patterns to review (e.g., script.py, *.py).")

    args = parser.parse_args()

    # If no arguments provided, show help
    if not args.files and not args.recursive:
        parser.print_help()
        sys.exit(1)

    files_to_review = get_files(args)

    if not files_to_review:
        print("No matching files found.")
        sys.exit(0)

    print(f"Found {len(files_to_review)} file(s) scheduled for review.")

    for i, f in enumerate(files_to_review):
        review_file(f)

if __name__ == "__main__":
    main()