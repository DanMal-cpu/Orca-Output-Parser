# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:10:16 2026

@author: dmalcolm
"""

import os

from contextlib import contextmanager
from glob import glob
from shutil import copy2
from subprocess import run




def cleanup_orca_files(orca_dir, base_filename):
    """ Remove .out and .gbw files from ORCA directory after processing. """
    
    for file in [f"{base_filename}.out", f"{base_filename}.gbw"]:
        file_path = os.path.join(orca_dir, file)                    # Searches for file in orca_dir.
        if os.path.exists(file_path):
            os.remove(file_path)



def cleanup_temp_files(orca_dir, base_filename):
    """ Optional: Clean up temporary files after organization. """
    
    temp_patterns = ["*.bub", "*.tmp"]                              # Patterns for temporary files to remove.
    removed_count = 0
    
    for pattern in temp_patterns:
        for file_path in glob(os.path.join(orca_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"   🗑️  Removed: {os.path.basename(file_path)}")
                removed_count += 1
            except FileNotFoundError:
                pass
            except PermissionError:
                print(f"   ⚠️  Permission denied: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"   ⚠️  Could not remove {os.path.basename(file_path)}: {e}")
    
    if removed_count > 0:
        print(f"\n   ✅ Removed {removed_count} temporary files")



def copy_orca_files_to_orca_dir(source_dir, target_dir, base_filename):
    """ Copy .out and .gbw files from source directory to ORCA directory. """
    
    files_to_copy = [f"{base_filename}.out", f"{base_filename}.gbw"]
    copied_files = []
    
    for file in files_to_copy:
        source_path = os.path.join(source_dir, file)                # python_directory.
        target_path = os.path.join(target_dir, file)                # orca_directory.
        
        if os.path.exists(source_path):
            copy2(source_path, target_path)                         # copy2 preserves metadata.
            print(f"   📋 Copied: {file} -> {target_dir}")
            copied_files.append(file)
    
    return len(copied_files) == len(files_to_copy)



def generate_molden_file(orca_dir, base_filename):
    """ Generate .molden.input file from ORCA .gbw file using orca_2mkl. """
    
    with working_directory(orca_dir):                                               # Check if .gbw file exists.
        gbw_file = f"{base_filename}.gbw"
        if not os.path.exists(gbw_file):
            print(f"❌ Error: {gbw_file} not found!")
            return False
        
        orca_2mkl_cmd = None                                                        # Check if orca_2mkl exists (it might be orca_2mkl.exe on Windows).
        possible_names = ["orca_2mkl.exe", "orca_2mkl", "orca_2mkl.bat"]
        
        for name in possible_names:
            if os.path.exists(name):
                orca_2mkl_cmd = name
                break
        
        if not orca_2mkl_cmd:
            print(f"❌ Error: orca_2mkl not found in {orca_dir}!")
            return False
        
        original_molden = f"{base_filename}.molden.input"                           # Remove existing molden files to ensure overwrite.
        new_molden = f"{base_filename}.molden"
        
        for file in [original_molden, new_molden]:    
            if os.path.exists(file):
                os.remove(file)
                print(f"\n   🔄 Removing existing {file}")
        
        cmd = [orca_2mkl_cmd, base_filename, "-molden"]
        
        try:
            result = run(cmd, capture_output=True, text=True, check=False)          # Run the command and capture output.
        
            if result.returncode == 0 and os.path.exists(original_molden):          # ORCA creates {filename}.molden.input.
                os.rename(original_molden, new_molden)                              # Change file type from .molden.input to .molden.
                size = os.path.getsize(new_molden)
                print(f"\n   ✅ {new_molden} generated ({size:,} bytes)\n")
                return True
            else:
                print("\n❌ Failed to generate molden file!")
                if result.stderr:
                    print(f"Error: {result.stderr.strip()}")
                return False
                
        except Exception as e:
            print(f"\n❌ Exception occurred: {str(e)}")
            return False
        


@contextmanager
def working_directory(path):
    """ Context manager for temporarily changing directory. """
    
    original_dir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(original_dir)