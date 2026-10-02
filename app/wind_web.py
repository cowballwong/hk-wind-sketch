"""HK Wind Sketch: the page as a Flask blueprint.

One code path, two hosts (same pattern as the HK landscape checker):
  * a personal dashboard can register it with its own base template;
  * app.py here registers it with a plain standalone base.

The page is static: the wind simulation runs in the visitor's browser (WebGL2,
or plain JavaScript), and the per-area data files are small JSON files in
static/data/. There is no API and nothing is uploaded.
"""
from __future__ import annotations

from pathlib import Path

from flask import Blueprint, render_template


def create_blueprint(*, base_template: str, ctx_fn=None) -> Blueprint:
    here = Path(__file__).resolve().parent
    bp = Blueprint("wind", __name__, template_folder=str(here / "templates"),
                   static_folder=str(here / "static"), static_url_path="/wind-static")

    @bp.route("/wind")
    def page():
        ctx = ctx_fn() if ctx_fn else {}
        return render_template("wind/page.html", base_template=base_template, **ctx)

    return bp
