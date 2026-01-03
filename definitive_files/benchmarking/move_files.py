import os
import shutil
import argparse

def move_files(to_move, destination_dir):
    # Ensure destination directory exists
    if not os.path.exists(destination_dir):
        print(f"Directory '{destination_dir}' does not exist. Creating it...")
        os.makedirs(destination_dir)

    print(f"\nFiles to move to '{destination_dir}':")
    # Filter list to only show files that actually exist to avoid confusion
    valid_files = [f for f in to_move if os.path.exists(f)]
    
    if not valid_files:
        print("  (No matching files found)")
        return

    for f in valid_files:
        print(" ", f)

    confirm = input("Move these files? [y/N]: ").lower()
    if confirm == "y":
        for f in valid_files:
            try:
                # shutil.move(source, destination)
                shutil.move(f, destination_dir)
                print(f"Moved: {f}")
            except Exception as e:
                print(f"Error moving {f}: {e}")
    else:
        print("Operation cancelled.")

# --- usage ---

parser = argparse.ArgumentParser(description="Script che recupera tutti i logs in questa cartella e realizza una tabella |process | times|")
parser.add_argument('-p', '--path', type=str, default="./scalability_folder", help='Percorso della directory da cui recuperare i logs (es. --path \\directory)')

args = parser.parse_args()
target_directory = args.path

# 1. Select files starting with prefixes
to_move_logs = [f for f in os.listdir(".") if f.startswith("benchmarking.sh.o") or f.startswith("benchmarking.sh.e")]
move_files(to_move_logs, target_directory)

# 2. Select specific filenames
to_move_results = ['results.txt', 'efficiency.png', 'speedup.png', 'table.txt']
move_files(to_move_results, target_directory)