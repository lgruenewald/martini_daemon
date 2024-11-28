#!/usr/bin/env python3

import os


def backup_try(path):
    if os.path.isfile(path):
        bkup_num = 2
        bkup_path = f"#{path}.1#"
        while os.path.isfile(bkup_path):
            bkup_num += 1
            bkup_path = f"#{path}.{bkup_num}#"
        os.rename(path, bkup_path)
        print(f"Backed up {path} to {bkup_path}")
