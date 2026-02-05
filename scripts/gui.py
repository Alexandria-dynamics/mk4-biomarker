#!/usr/bin/env python3
"""
MK4 Biomarker - GUI
Alexandria Dynamics

Simple GUI wrapper around parse_geo_soft.py.
Select a .soft.gz file, run analysis, open HTML report.
"""

import shutil
import threading
import webbrowser
from pathlib import Path
from tkinter import Tk, Label, Button, filedialog, StringVar, Frame

# Import the analysis pipeline
from parse_geo_soft import (
    INPUT_DIR,
    parse_soft_file,
    compute_chaos_metrics,
    create_output_directory,
    generate_html_report,
    PROJECT_ROOT,
)


class MK4App:
    def __init__(self, root):
        self.root = root
        self.root.title("MK4 Biomarker Analysis")
        self.root.geometry("500x300")
        self.root.resizable(False, False)
        self.root.configure(bg="#f5f7fa")

        self.status = StringVar(value="Select a .soft.gz file to begin.")
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

        self.btn_select = Button(
            content,
            text="Select .soft.gz file",
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
            text="Open Report",
            font=("Arial", 13),
            bg="#2ecc71",
            fg="white",
            activebackground="#27ae60",
            activeforeground="white",
            padx=20,
            pady=10,
            relief="flat",
            cursor="hand2",
            command=self.open_report,
        )
        # Hidden until report is ready
        self.btn_open.pack(pady=(0, 10))
        self.btn_open.pack_forget()

    def select_file(self):
        filepath = filedialog.askopenfilename(
            title="Select GEO SOFT file",
            filetypes=[
                ("GEO SOFT files", "*.soft.gz *.soft"),
                ("All files", "*.*"),
            ],
        )
        if not filepath:
            return

        filepath = Path(filepath)
        self.status.set(f"Selected: {filepath.name}")
        self.btn_open.pack_forget()
        self.btn_select.configure(state="disabled")

        # Run analysis in background thread to keep GUI responsive
        thread = threading.Thread(target=self.run_analysis, args=(filepath,), daemon=True)
        thread.start()

    def run_analysis(self, filepath):
        try:
            # Copy file to data/input/ so paths stay consistent
            dest = INPUT_DIR / filepath.name
            if dest != filepath:
                INPUT_DIR.mkdir(parents=True, exist_ok=True)
                shutil.copy2(filepath, dest)

            self.status.set("Parsing file...")
            dataset = parse_soft_file(dest)

            if len(dataset["samples"]) == 0:
                self.status.set("No valid samples found in this file.")
                self.btn_select.configure(state="normal")
                return

            self.status.set("Computing chaos metrics...")
            chaos_results = compute_chaos_metrics(dataset)

            self.status.set("Generating HTML report...")
            output_dir = create_output_directory(dataset["dataset_name"])
            html_path = generate_html_report(output_dir, dataset, chaos_results)

            self.report_path = html_path
            self.status.set(f"Done! Report: {html_path.relative_to(PROJECT_ROOT)}")
            self.btn_open.pack(pady=(0, 10))

        except Exception as e:
            self.status.set(f"Error: {e}")

        finally:
            self.btn_select.configure(state="normal")

    def open_report(self):
        if self.report_path and self.report_path.exists():
            webbrowser.open(self.report_path.as_uri())


def main():
    root = Tk()
    MK4App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
