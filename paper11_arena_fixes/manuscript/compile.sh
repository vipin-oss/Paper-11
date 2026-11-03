#!/usr/bin/env bash
# =============================================================================
# Paper 11 - one-touch LaTeX build script
#
# REQUIRED TOOLS (both must be on PATH):
#     pdflatex   - PDFLaTeX   (texlive-latex-base / texlive-latex-recommended)
#     bibtex     - BibTeX     (texlive-binaries / texlive-latex-base)
#
# The script fails with a clear message if either tool is missing. It never
# silently reuses a previously generated main.pdf: the PDF in this directory is
# always the output of THIS build.
#
# Debian / Ubuntu:
#     sudo apt-get install texlive-latex-base texlive-latex-recommended \
#                          texlive-fonts-recommended texlive-latex-extra
# =============================================================================
set -euo pipefail

echo "=== Paper 11: clean LaTeX build (pdflatex + bibtex) ==="

# ---------------------------------------------------------------- dependencies
command -v pdflatex >/dev/null 2>&1 || {
    echo "ERROR: pdflatex is required."
    echo "Install a TeX distribution, e.g.:"
    echo "    sudo apt-get install texlive-latex-base texlive-latex-recommended \\"
    echo "                         texlive-fonts-recommended texlive-latex-extra"
    echo "No PDF was produced or modified."
    exit 1
}

command -v bibtex >/dev/null 2>&1 || {
    echo "ERROR: BibTeX is required for a clean bibliography rebuild."
    echo "Install a TeX distribution, e.g.:"
    echo "    sudo apt-get install texlive-latex-base texlive-latex-recommended \\"
    echo "                         texlive-fonts-recommended texlive-latex-extra"
    echo "No PDF was produced or modified."
    exit 1
}

echo "  pdflatex : $(command -v pdflatex)"
echo "  bibtex   : $(command -v bibtex)"

# ---------------------------------------------------------------- clean build
rm -f main.aux main.log main.out main.blg main.spl main.toc

pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex

# ---------------------------------------------------------------- report
if [ ! -f main.pdf ]; then
    echo "ERROR: pdflatex finished but main.pdf was not created." >&2
    exit 1
fi

UNDEF=$(grep -c "LaTeX Warning: Reference" main.log 2>/dev/null || true)
OVER=$(grep -c "Overfull" main.log 2>/dev/null || true)

echo
echo "=== Compilation complete ==="
ls -lh main.pdf
echo "  undefined references : ${UNDEF:-0}"
echo "  overfull boxes       : ${OVER:-0}"
if [ "${UNDEF:-0}" -ne 0 ] || [ "${OVER:-0}" -ne 0 ]; then
    echo "  (warnings present - inspect main.log; these do not stop the build)"
fi
echo "=== main.pdf generated from main.tex by this run ==="
