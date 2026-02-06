#!/usr/bin/env python3
"""
MK4 Biomarker - GUI (OPRAVENÝ)
Alexandria Dynamics

Simple GUI wrapper around parse_geo_soft.py.
Supports BOTH SOFT and MATRIX formats!

Changes:
- Uses parse_file() instead of parse_soft_file()
- Supports .txt.gz and .txt files (matrix format)
- Auto-detects format
"""

import shutil
import threading
import webbrowser
from pathlib import Path
from tkinter import Tk, Label, Button, filedialog, StringVar, Frame

# Import the analysis pipeline (UPDATED!)
from parse_geo_soft import (
    INPUT_DIR,
    parse_file,  # ← ZMĚNA: unified parser místo parse_soft_file
    compute_chaos_metrics,
    create_output_directory,
    save_results,
    generate_plot,
    PROJECT_ROOT,
)


class MK4App:
    def __init__(self, root):
        self.root = root
        self.root.title("MK4 Biomarker Analysis")
        self.root.geometry("500x350")
        self.root.resizable(False, False)
        self.root.configure(bg="#f5f7fa")

        self.status = StringVar(value="Select a file to begin.")
        self.report_path = None

        # Header
        header = Frame(root, bg="#667eea", height=60)
        header.pack(fill="x")
        header.pack_propagate(False)
        Label(
            header,
            text="MK4 Biomarker Analysis",
            font=("Arial", 18, "bold"),
            bg="#667eea",
            fg="white",
        ).pack(expand=True)

        # Content
        content = Frame(root, bg="#f5f7fa", padx=30, pady=20)
        content.pack(fill="both", expand=True)

        # Info label
        Label(
            content,
            text="Supports: SOFT format (.soft.gz, .soft)\nand MATRIX format (.txt.gz, .txt)",
            font=("Arial", 9),
            bg="#f5f7fa",
            fg="#888",
        ).pack(pady=(5, 10))

        self.btn_select = Button(
            content,
            text="Select File",
            font=("Arial", 13),
            bg="#667eea",
            fg="white",
            activebackground="#764ba2",
            activeforeground="white",
            padx=20,
            pady=10,
            relief="flat",
            cursor="hand2",
            command=self.select_file,
        )
        self.btn_select.pack(pady=(10, 15))

        self.lbl_status = Label(
            content,
            textvariable=self.status,
            font=("Arial", 11),
            bg="#f5f7fa",
            fg="#555",
            wraplength=440,
        )
        self.lbl_status.pack(pady=(0, 15))

        self.btn_open = Button(
            content,
            text="Open Results Folder",
            font=("Arial", 13),
            bg="#2ecc71",
            fg="white",
            activebackground="#27ae60",
            activeforeground="white",
            padx=20,
            pady=10,
            relief="flat",
            cursor="hand2",
            command=self.open_results,
        )
        # Hidden until analysis is done
        self.btn_open.pack(pady=(0, 10))
        self.btn_open.pack_forget()

    def select_file(self):
        filepath = filedialog.askopenfilename(
            title="Select GEO data file",
            filetypes=[
                ("GEO SOFT files", "*.soft.gz *.soft"),
                ("Matrix files", "*.txt.gz *.txt"),  # ← PŘIDÁNO!
                ("All files", "*.*"),
            ],
        )
        if not filepath:
            return

        filepath = Path(filepath)
        self.status.set(f"Selected: {filepath.name}")
        self.btn_open.pack_forget()
        self.btn_select.configure(state="disabled")

        # Run analysis in background thread
        thread = threading.Thread(target=self.run_analysis, args=(filepath,), daemon=True)
        thread.start()

    def run_analysis(self, filepath):
        try:
            # Copy file to data/input/
            dest = INPUT_DIR / filepath.name
            if dest != filepath:
                INPUT_DIR.mkdir(parents=True, exist_ok=True)
                shutil.copy2(filepath, dest)

            self.status.set("Detecting file format...")
            
            # ═══════════════════════════════════════════════════════
            # ZMĚNA: Použij parse_file() místo parse_soft_file()
            # ═══════════════════════════════════════════════════════
            dataset = parse_file(dest)
            
            if dataset is None:
                self.status.set("Unknown file format. Supported: SOFT (.soft.gz) or MATRIX (.txt.gz)")
                self.btn_select.configure(state="normal")
                return

            if len(dataset["samples"]) == 0:
                self.status.set("No valid samples found in this file.")
                self.btn_select.configure(state="normal")
                return

            self.status.set("Computing chaos metrics...")
            chaos_results = compute_chaos_metrics(dataset)

            self.status.set("Saving results...")
            output_dir = create_output_directory(dataset["dataset_name"])
            save_results(output_dir, dataset, chaos_results)
            generate_plot(output_dir, chaos_results)

            self.report_path = output_dir
            
            # Summary
            if 'error' not in chaos_results['statistics']:
                stats = chaos_results['statistics']
                summary = (
                    f"✓ Analysis complete!\n"
                    f"Control: {stats['control_mean']:.3f} ± {stats['control_std']:.3f}\n"
                    f"Disease: {stats['disease_mean']:.3f} ± {stats['disease_std']:.3f}\n"
                    f"p-value: {stats['p_value']:.2e} "
                    f"{'***' if stats['p_value'] < 0.001 else '**' if stats['p_value'] < 0.01 else '*' if stats['p_value'] < 0.05 else 'ns'}"
                )
            else:
                summary = "✓ Analysis complete! (insufficient samples for statistics)"
            
            self.status.set(summary)
            self.btn_open.pack(pady=(0, 10))

        except Exception as e:
            self.status.set(f"Error: {e}")
            import traceback
            traceback.print_exc()

        finally:
            self.btn_select.configure(state="normal")

    def open_results(self):
        if self.report_path and self.report_path.exists():
            # Open folder in file manager
            import os
            import platform
            
            system = platform.system()
            if system == 'Darwin':  # macOS
                os.system(f'open "{self.report_path}"')
            elif system == 'Windows':
                os.startfile(self.report_path)
            else:  # Linux
                os.system(f'xdg-open "{self.report_path}"')


def main():
    root = Tk()
    app = MK4App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
