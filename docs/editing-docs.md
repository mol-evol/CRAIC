# Editing these docs

The documentation is plain **Markdown** in `docs/`, rendered by **MkDocs** with the
**Material** theme. There is no database and no proprietary format to learn — edit a file,
refresh, done.

## Live preview

```bash
./serve-docs.command        # installs mkdocs-material into .venv, then serves on http://127.0.0.1:8000
# or, if you prefer:
.venv/bin/python -m mkdocs serve
```

The site rebuilds automatically on every save.

## Add or change a page

1. Create a Markdown file under `docs/` (for example `docs/guide/my-topic.md`).
2. Add it to the `nav:` tree in `mkdocs.yml`.
3. Write Markdown. Admonitions, content tabs, tables, footnotes, and syntax highlighting are
   all enabled.

!!! example "Callouts"
    Lines beginning `!!! note`, `!!! tip`, `!!! warning`, `!!! example` render as boxes like
    this one.

## Images and video

Put assets in `docs/assets/` and reference them relative to the page:

```markdown
![A screenshot](assets/my-image.png){ width="100%" }
```

For local video, drop an `.mp4` in `docs/assets/` and embed it with HTML (Markdown allows raw
HTML):

```html
<video controls width="100%" poster="assets/poster.png">
  <source src="assets/walkthrough.mp4" type="video/mp4">
</video>
```

External video works too — paste a YouTube/Vimeo `<iframe>` embed.

## Build a static site

```bash
.venv/bin/python -m mkdocs build     # produces a self-contained site/ folder
```

`site/` is plain HTML/CSS/JS: open it in any browser, or host it anywhere — GitHub Pages,
Read the Docs, an S3 bucket, or an internal server. Nothing about the output is tied to a
platform.

## Updating the citation

While the paper is under review, CRAIC is cited as the bioRxiv preprint. The citation
is written out in five places; change all of them together (the test suite checks that
the first three agree):

1. `craic/__init__.py`: `CITATION`, which the app's About box and copy button use;
2. `README.md`: the **Citation** section;
3. `docs/about.md`: the citation and the BibTeX entry below it;
4. `docs/index.md`: the **Citing CRAIC** box;
5. `CITATION.cff`: the `preferred-citation` block, which drives GitHub's **Cite this
   repository** button.

**When bioRxiv posts the preprint**, replace "link to follow" with its DOI link
(`https://doi.org/…`), and add `doi: "…"` to `preferred-citation` in `CITATION.cff`.

**When the paper is published**, replace the preprint with the journal citation
(journal, volume, pages, DOI), and in `CITATION.cff` set `journal`, `volume`, `start`,
`end` and `doi` to match.

Push to `main` and the website updates itself.

## Why this stack

Markdown keeps the *source* trivially editable and diff-friendly; Material renders it into a
rich, searchable, themed site; the build is a static folder, so hosting is unconstrained. See
the handover notes for the full rationale, or [Architecture & development](development.md) for
how to turn on auto-generated API docs from the code's docstrings.
