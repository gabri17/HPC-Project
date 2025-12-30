import os
import re

def delete_files(to_delete):

    print("Files to delete:")
    for f in to_delete:
        print(" ", f)

    confirm = input("Delete these files? [y/N]: ").lower()
    if confirm == "y":
        for f in to_delete:
            try:
                os.remove(f)
                print("Deleted:", f)
            except FileNotFoundError:
                print("Not exists:", f) 

to_delete = [f for f in os.listdir(".") if f.startswith("benchmarking.sh.o") or f.startswith("benchmarking.sh.e")]
delete_files(to_delete)
to_delete = ['results.txt', 'efficiency.png', 'speedup.png', 'table.txt']
delete_files(to_delete)