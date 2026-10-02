"""HK Wind Sketch: standalone host.

    pip install flask
    python app.py        then open http://127.0.0.1:5058/wind

Or skip Python entirely: open hk-wind-sketch.html (everything inlined).
"""
from __future__ import annotations

import os

from flask import Flask, redirect

import wind_web

app = Flask(__name__)
app.register_blueprint(wind_web.create_blueprint(base_template="wind/standalone_base.html"))


@app.route("/")
def home():
    return redirect("/wind")


if __name__ == "__main__":
    app.run(host=os.environ.get("HKW_HOST", "127.0.0.1"), port=int(os.environ.get("HKW_PORT", "5058")), debug=False)
