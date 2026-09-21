"""The one HTML page both reports are: `kmerseek_hits_template.html` with its two tokens
filled. The search report passes the whole result; the pair figure passes `{"pair":
model}`, which the template renders as a one-row page with the row open. Everything on
the page is drawn by the template's own script from the JSON, so the Python side never
writes markup."""

import html
import json
import os

TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kmerseek_hits_template.html")


def embed_json(value):
    """JSON safe inside a <script>: `</` cannot close the tag."""
    return json.dumps(value).replace("</", "<\\/")


def render_page(title, data):
    with open(TEMPLATE) as fh:
        page = fh.read()
    assert page.count("__TITLE__") == 1 and page.count("__DATA__") == 1
    return page.replace("__TITLE__", html.escape(title)).replace("__DATA__", embed_json(data))
