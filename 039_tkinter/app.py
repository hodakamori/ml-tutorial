import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
import pandas as pd
from sklearn.linear_model import LinearRegression

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class RegressionApp:
    def __init__(self, master):
        self.master = master
        master.title("線形回帰アプリ")

        self.frame = ttk.Frame(master, padding=10)
        self.frame.pack(fill=tk.BOTH, expand=True)

        self.file_label = ttk.Label(self.frame, text="CSVファイル:")
        self.file_label.grid(row=0, column=0, sticky=tk.E, pady=5)

        self.file_path = tk.StringVar()
        self.file_entry = ttk.Entry(self.frame, textvariable=self.file_path, width=50)
        self.file_entry.grid(row=0, column=1, pady=5)

        self.browse_button = ttk.Button(
            self.frame, text="ファイルを選択", command=self.browse_file
        )
        self.browse_button.grid(row=0, column=2, padx=5, pady=5)

        self.run_button = ttk.Button(
            self.frame, text="実行", command=self.run_regression
        )
        self.run_button.grid(row=1, column=2, padx=5, pady=5)

        self.canvas = None

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="CSVファイルを選択",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if filename:
            self.file_path.set(filename)

    def run_regression(self):
        file = self.file_path.get()
        if not file:
            messagebox.showerror("エラー", "CSVファイルを選択してください。")
            return

        try:
            df = pd.read_csv(file)
        except Exception as e:
            messagebox.showerror("エラー", f"ファイルの読み込みに失敗しました:\n{e}")
            return

        if df.shape[1] < 2:
            messagebox.showerror("エラー", "データが2列以上必要です( x と y )。")
            return

        X = df.iloc[:, 0].values.reshape(-1, 1)
        y = df.iloc[:, 1].values

        model = LinearRegression()
        model.fit(X, y)

        y_pred = model.predict(X)

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(X, y, color="blue", label="Data")
        ax.plot(X, y_pred, color="red", linewidth=2, label="Regression Line")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_title("Linear Regression Result")
        ax.legend()

        if self.canvas:
            self.canvas.get_tk_widget().destroy()

        self.canvas = FigureCanvasTkAgg(fig, master=self.frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=2, column=0, columnspan=3, pady=10)


if __name__ == "__main__":
    root = tk.Tk()
    app = RegressionApp(root)
    root.mainloop()
