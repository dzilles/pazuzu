#!/usr/bin/env python3
import argparse
import glob
import os
import shutil
import subprocess
import sys

# --- Configuration ---
GEMINI_CMD = "gemini"

# The prompt for deriving unit tests
TEST_DERIVATION_PROMPT = """
You are an expert Software Quality Assurance Engineer and Python Developer.
Your task is to analyze the provided source code file and derive a comprehensive unit test plan.

**DO NOT IMPLEMENT THE TESTS YET.**
Instead, generate detailed documentation for the tests that *should* be implemented.

**Requirements:**
1.  **Coverage**: EVERY function and method in the source file must have at least one corresponding unit test case.
2.  **Naming**:
    *   Test Function Name: `test_<original_function_name>_<scenario>`
    *   Link each test explicitly to the function it targets.
3.  **Structure per Test Case**:
    *   **Target Function**: The name of the function being tested.
    *   **Test Name**: The proposed name for the test function.
    *   **Description**: What aspect/logic is being tested?
    *   **Input/Setup**: What arguments or state setup are required?
    *   **Expected Output/Behavior**: What should the return value or state change be?
    *   **Edge Cases**: Mention any specific edge cases (nulls, boundaries, errors) covered.

**Output Format:**
Please provide the output in Markdown format.
"""

def get_files(args):
    """Collects a list of files based on arguments (recursive or glob patterns)."""
    files = []
    if args.recursive:
        for root, _, filenames in os.walk("."):
            if any(x in root for x in [".git", "__pycache__", ".venv", ".vscode", "output"]):
                continue
            for f in filenames:
                files.append(os.path.join(root, f))
    elif args.files:
        for pattern in args.files:
            matched = glob.glob(pattern)
            if matched:
                files.extend(matched)
            else:
                files.append(pattern)
    return sorted(list(set(files)))

def derive_tests_for_file(file_path):
    """Starts a Gemini CLI instance to derive tests for the specified file."""
    if not os.path.isfile(file_path):
        print(f"Skipping '{file_path}': Not a file.")
        return

    print(f"\n{'='*60}")
    print(f"DERIVING TESTS FOR: {file_path}")
    print(f"{ '='*60}")

    full_prompt = f"{TEST_DERIVATION_PROMPT}\n\nPlease analyze the file: {file_path}\nIMPORTANT: Provide the test plan as text output only. Do NOT attempt to create or modify any files."

    try:
        # Resolve the executable path (handles .cmd/.bat on Windows)
        executable = shutil.which(GEMINI_CMD) or GEMINI_CMD
        
        print(f"Launching {executable} for {file_path}...")
        # Run in one-shot mode, capturing output.
        result = subprocess.run([executable, full_prompt], capture_output=True, text=True)
        
        if result.returncode == 0:
            # Construct path inside the output/ directory
            rel_path = os.path.normpath(file_path).lstrip(os.sep)
            output_filename = os.path.join("output", f"{rel_path}.tests.md")
            
            # Ensure the target subdirectory exists
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)
            
            with open(output_filename, "w") as f:
                f.write(result.stdout)
            print(f"Test plan saved to: {output_filename}")
        else:
            print(f"Error processing {file_path}:")
            print(result.stderr)
        
    except FileNotFoundError:
        print(f"Error: Command '{GEMINI_CMD}' not found. Please ensure it is installed and in your PATH.")
    except KeyboardInterrupt:
        print("\nSession interrupted.")

def main():
    parser = argparse.ArgumentParser(description="Derive unit test plans for files using the Gemini CLI.")
    parser.add_argument("-r", "--recursive", action="store_true", help="Recursively process all files in the current directory.")
    parser.add_argument("files", nargs="*", help="Files or glob patterns to process (e.g., script.py, *.py).")

    args = parser.parse_args()

    if not args.files and not args.recursive:
        parser.print_help()
        sys.exit(1)

    files_to_process = get_files(args)

    if not files_to_process:
        print("No matching files found.")
        sys.exit(0)

    print(f"Found {len(files_to_process)} file(s) scheduled for test derivation.")

    for i, f in enumerate(files_to_process):
        derive_tests_for_file(f)

if __name__ == "__main__":
    main()
