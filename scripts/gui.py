#!/usr/bin/env python3
"""
MK4 Biomarker - GUI
Alexandria Dynamics

Simple GUI for batch processing of GEO files.
Features:
- Select multiple files at once
- Each dataset gets its own subdirectory
- Error logging to error.log
- Completed files moved to complete/ folder
- Failed files moved to error/ folder
- Processing continues even if one file fails
"""

import shutil
import threading
import traceback
import webbrowser
from datetime import datetime
from pathlib import Path
from tkinter import Tk, Label, Button, filedialog, StringVar, Frame, Listbox, Scrollbar, END, VERTICAL, BOTH, LEFT, RIGHT, Y

# Import from auto_parse_geo (universal parser)
from auto_parse_geo import (
    INPUT_DIR,
    RESULTS_DIR,
    PROJECT_ROOT,
    auto_analyze,
    detect_file_format,
)


class MK4BatchApp:
    def __init__(self, root):
        self.root = root
        self.root.title("MK4 Biomarker Analysis")
        self.root.geometry("700x500")
        self.root.configure(bg="#f5f7fa")

        self.selected_files = []
        self.run_dir = None
        self.complete_dir = None
        self.error_dir = None
        self.log_file = None

        self.status = StringVar(value="Select files to begin batch processing.")

        # Header
        header = Frame(root, bg="#667eea", height=50)
        header.pack(fill="x")
        header.pack_propagate(False)
        Label(
            header,
            text="MK4 Biomarker - Batch Analysis",
            font=("Arial", 16, "bold"),
            bg="#667eea",
            fg="white",
        ).pack(expand=True)

        # Main content
        content = Frame(root, bg="#f5f7fa", padx=20, pady=15)
        content.pack(fill=BOTH, expand=True)

        # Button row
        btn_frame = Frame(content, bg="#f5f7fa")
        btn_frame.pack(fill="x", pady=(0, 10))

        self.btn_select = Button(
            btn_frame,
            text="Select Files",
            font=("Arial", 11),
            bg="#667eea",
            fg="white",
            activebackground="#764ba2",
            activeforeground="white",
            padx=15,
            pady=5,
            relief="flat",
            cursor="hand2",
            command=self.select_files,
        )
        self.btn_select.pack(side=LEFT, padx=(0, 10))

        self.btn_run = Button(
            btn_frame,
            text="Run Analysis",
            font=("Arial", 11),
            bg="#2ecc71",
            fg="white",
            activebackground="#27ae60",
            activeforeground="white",
            padx=15,
            pady=5,
            relief="flat",
            cursor="hand2",
            command=self.run_analysis,
            state="disabled",
        )
        self.btn_run.pack(side=LEFT, padx=(0, 10))

        self.btn_open = Button(
            btn_frame,
            text="Open Results",
            font=("Arial", 11),
            bg="#3498db",
            fg="white",
            activebackground="#2980b9",
            activeforeground="white",
            padx=15,
            pady=5,
            relief="flat",
            cursor="hand2",
            command=self.open_results,
            state="disabled",
        )
        self.btn_open.pack(side=LEFT)

        # File list with scrollbar
        list_frame = Frame(content, bg="#f5f7fa")
        list_frame.pack(fill=BOTH, expand=True, pady=(0, 10))

        scrollbar = Scrollbar(list_frame, orient=VERTICAL)
        self.file_list = Listbox(
            list_frame,
            font=("Consolas", 10),
            selectmode="extended",
            yscrollcommand=scrollbar.set,
            bg="white",
            relief="flat",
            highlightthickness=1,
            highlightbackground="#ddd",
        )
        scrollbar.config(command=self.file_list.yview)
        scrollbar.pack(side=RIGHT, fill=Y)
        self.file_list.pack(side=LEFT, fill=BOTH, expand=True)

        # Status label
        self.lbl_status = Label(
            content,
            textvariable=self.status,
            font=("Arial", 10),
            bg="#f5f7fa",
            fg="#555",
            anchor="w",
            wraplength=650,
        )
        self.lbl_status.pack(fill="x", pady=(5, 0))

        # Progress label
        self.progress_var = StringVar(value="")
        self.lbl_progress = Label(
            content,
            textvariable=self.progress_var,
            font=("Arial", 10, "bold"),
            bg="#f5f7fa",
            fg="#667eea",
            anchor="w",
        )
        self.lbl_progress.pack(fill="x")

    def select_files(self):
        """Open file dialog to select multiple files."""
        filepaths = filedialog.askopenfilenames(
            title="Select GEO files",
            filetypes=[
                ("GEO files", "*.soft.gz *.soft *.txt.gz *.txt"),
                ("SOFT files", "*.soft.gz *.soft"),
                ("Matrix files", "*.txt.gz *.txt"),
                ("All files", "*.*"),
            ],
        )

        if filepaths:
            self.selected_files = [Path(f) for f in filepaths]
            self.file_list.delete(0, END)
            for f in self.selected_files:
                self.file_list.insert(END, f.name)

            self.status.set(f"Selected {len(self.selected_files)} file(s). Click 'Run Analysis' to start.")
            self.btn_run.config(state="normal")
            self.btn_open.config(state="disabled")

    def run_analysis(self):
        """Start batch analysis in background thread."""
        if not self.selected_files:
            return

        self.btn_select.config(state="disabled")
        self.btn_run.config(state="disabled")
        self.status.set("Starting batch analysis...")

        thread = threading.Thread(target=self._process_batch, daemon=True)
        thread.start()

    def _process_batch(self):
        """Process all selected files (runs in background thread)."""
        import re
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Create main run directory
        self.run_dir = RESULTS_DIR / f"batch_{timestamp}"
        self.run_dir.mkdir(parents=True, exist_ok=True)

        # Create complete and error directories in INPUT_DIR
        self.complete_dir = INPUT_DIR / "complete"
        self.error_dir = INPUT_DIR / "error"
        self.complete_dir.mkdir(exist_ok=True)
        self.error_dir.mkdir(exist_ok=True)

        # Create error log
        self.log_file = self.run_dir / "error.log"

        # Group files by GSE ID to combine SOFT + Matrix pairs
        file_groups = {}
        for filepath in self.selected_files:
            # Extract GSE ID from filename
            match = re.search(r'(GSE\d+)', filepath.name, re.IGNORECASE)
            gse_id = match.group(1).upper() if match else filepath.stem

            if gse_id not in file_groups:
                file_groups[gse_id] = {"soft": None, "matrix": None, "files": []}

            file_groups[gse_id]["files"].append(filepath)

            # Detect file type
            fmt = detect_file_format(filepath)
            if fmt == "SOFT":
                file_groups[gse_id]["soft"] = filepath
            elif fmt == "MATRIX":
                file_groups[gse_id]["matrix"] = filepath

        total = len(file_groups)
        completed = 0
        failed = 0
        file_idx = 0

        for gse_id, group in file_groups.items():
            files = group["files"]
            self.progress_var.set(f"Processing {gse_id} ({len(files)} file(s))")

            try:
                # Process as a group (auto_analyze handles combining)
                self._process_file_group(group)
                completed += 1

                # Move all files in group to complete
                for filepath in files:
                    try:
                        dest = self.complete_dir / filepath.name
                        if dest.exists():
                            dest = self.complete_dir / f"{filepath.stem}_{timestamp}{filepath.suffix}"
                        shutil.move(str(filepath), str(dest))
                    except Exception as move_err:
                        self._log_error(filepath.name, f"Move error: {move_err}")

                # Update list items
                for filepath in files:
                    idx = self.selected_files.index(filepath)
                    self._update_list_item(idx, f"[OK] {filepath.name}")

            except Exception as e:
                failed += 1
                error_msg = f"{type(e).__name__}: {e}"
                for filepath in files:
                    self._log_error(filepath.name, error_msg)
                self._log_error(gse_id, traceback.format_exc())

                # Move all files in group to error
                for filepath in files:
                    try:
                        dest = self.error_dir / filepath.name
                        if dest.exists():
                            dest = self.error_dir / f"{filepath.stem}_{timestamp}{filepath.suffix}"
                        shutil.move(str(filepath), str(dest))
                    except Exception as move_err:
                        self._log_error(filepath.name, f"Move error: {move_err}")

                    idx = self.selected_files.index(filepath)
                    self._update_list_item(idx, f"[ERROR] {filepath.name}")

                continue

        # Done
        self.progress_var.set("")
        self.status.set(
            f"Batch complete: {completed} succeeded, {failed} failed. "
            f"Results in: {self.run_dir.name}"
        )
        self.btn_select.config(state="normal")
        self.btn_open.config(state="normal")

    def _process_file_group(self, group):
        """Process a group of related files (SOFT + Matrix with same GSE ID)."""
        files = group["files"]
        soft_file = group.get("soft")
        matrix_file = group.get("matrix")

        # Copy files to input dir
        work_files = []
        for filepath in files:
            if filepath.parent != INPUT_DIR:
                dest = INPUT_DIR / filepath.name
                shutil.copy2(filepath, dest)
                work_files.append(dest)
            else:
                work_files.append(filepath)

        # Call auto_analyze with all files in the group
        result = auto_analyze(work_files)

        # Check for errors
        if "error" in result:
            raise ValueError(result["error"])

        dataset = result.get("dataset", {})
        analysis = result.get("results", {})

        # Check if we have expression data
        samples = dataset.get("samples", {})
        n_with_expr = sum(1 for s in samples.values() if s.get("expression"))
        if n_with_expr == 0:
            file_names = ", ".join(f.name for f in files)
            raise ValueError(f"No expression data found in: {file_names}")

        if "error" in analysis and analysis.get("n_control", 0) == 0:
            raise ValueError(f"No labeled samples found")

        # Move output to our batch directory structure
        output_dir_str = result.get("output_dir", "")
        if output_dir_str and output_dir_str != ".":
            auto_output = Path(output_dir_str)
            if auto_output.exists() and auto_output.is_dir() and auto_output != Path("."):
                # Move contents to our batch subdirectory
                dataset_name = dataset.get("dataset_name", files[0].stem)
                dataset_dir = self.run_dir / dataset_name

                if dataset_dir.exists():
                    shutil.rmtree(dataset_dir)

                shutil.move(str(auto_output), str(dataset_dir))

    def _log_error(self, filename, message):
        """Append error to log file."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self.log_file, "a", encoding="utf-8") as f:
            f.write(f"[{timestamp}] {filename}\n{message}\n\n")

    def _update_list_item(self, index, text):
        """Update listbox item (thread-safe via root.after)."""
        def update():
            self.file_list.delete(index)
            self.file_list.insert(index, text)
        self.root.after(0, update)

    def open_results(self):
        """Open results directory or index.html in browser."""
        if self.run_dir and self.run_dir.exists():
            # Try to find an index.html in any subdirectory
            for subdir in self.run_dir.iterdir():
                if subdir.is_dir():
                    index_file = subdir / "index.html"
                    if index_file.exists():
                        webbrowser.open(index_file.as_uri())
                        return

            # Fallback: open the directory
            webbrowser.open(self.run_dir.as_uri())


def main():
    root = Tk()
    MK4BatchApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
