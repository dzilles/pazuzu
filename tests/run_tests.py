import os
import subprocess
import sys

def run_tests():
    root_tests_dir = os.path.dirname(__file__)
    verification_dir = os.path.join(root_tests_dir, "verification")
    
    test_scripts = []
    for root, dirs, files in os.walk(verification_dir):
        if "verify.py" in files:
            test_scripts.append(os.path.join(root, "verify.py"))
    
    if not test_scripts:
        print("No tests found.")
        return

    print(f"Found {len(test_scripts)} tests.")
    
    failed_tests = []
    
    for script in test_scripts:
        print(f"\n--- Running {script} ---")
        try:
            # We don't capture output here so we can see progress
            subprocess.run([sys.executable, script], check=True)
        except subprocess.CalledProcessError:
            failed_tests.append(script)
            print(f"!!! Test {script} FAILED !!!")
    
    print("\n==============================")
    if failed_tests:
        print(f"Tests failed: {len(failed_tests)}")
        for t in failed_tests:
            print(f" - {t}")
        sys.exit(1)
    else:
        print("All tests PASSED.")
        sys.exit(0)

if __name__ == "__main__":
    run_tests()
