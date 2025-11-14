# ============================================================================
# JAVASCRIPT SEQUENCE HIGHLIGHTING FUNCTIONS
# ============================================================================

"""
JavaScript functions for sequence highlighting modal.

Mustache variables:
- `sequence_file_paths`: Array of FASTA file paths

Functions:
- `openHighlightPage()`: Open sequence highlighting modal
- `parseFasta()`: Parse FASTA format files
- `parseCsv()`: Parse CSV position data
- `generateHtml()`: Generate highlighted sequence HTML
- `highlightSequence()`: Apply color highlighting to sequences
- `loadFile()`: Fetch file content via HTTP
"""
script_sequence_str = raw"""
// --------------------------------------- sequence hightlighting part  ---------------------------------------------

let scrollPosition = 0;

function loadFile(filePath) {
    return fetch(filePath).then(response => response.text());
}

function openHighlightPage(csvFile) {
    // Store the current scroll position
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;

    const fastaFiles = [
        {{#:sequence_file_paths}}
        '{{{.}}}'{{^.[end]}}, {{/.[end]}}
        {{/:sequence_file_paths}}
    ];

    // Load all FASTA files and the CSV file
    const fastaPromises = fastaFiles.map(file => loadFile(file));
    Promise.all([...fastaPromises, loadFile(csvFile)])
        .then(contents => {
            const fastaContents = contents.slice(0, fastaFiles.length);
            const csvContent = contents[contents.length - 1];

            // Combine sequences from all FASTA files
            const sequences = fastaContents.flatMap(fastaContent => parseFasta(fastaContent));
            const highlights = parseCsv(csvContent);

            const htmlContent = generateHtml(sequences, highlights);

            // Insert the generated HTML into the modal
            document.getElementById('highlightedSequences').innerHTML = htmlContent;

            // Show the modal
            document.getElementById('highlightModal').style.display = "block";
        });
}

function parseFasta(fastaContent) {
    const lines = fastaContent.split('\n');
    const sequences = [];
    let currentHeader = '';
    let currentSequence = '';
    lines.forEach(line => {
        if (line.startsWith('>')) {
            if (currentHeader) {
                sequences.push({ header: currentHeader, sequence: currentSequence });
                currentSequence = '';
            }
            currentHeader = line;
        } else {
            currentSequence += line.trim();
        }
    });
    if (currentHeader) {
        sequences.push({ header: currentHeader, sequence: currentSequence });
    }
    return sequences;
}

function parseCsv(csvContent) {
    const lines = csvContent.trim().split('\n');
    const highlights = [];
    lines.slice(1).forEach(line => {
        const [seqIndex, startPosition, endPosition, iscomp] = line.split(',').map(Number);
        if (!isNaN(seqIndex) && !isNaN(startPosition) && !isNaN(endPosition) && !isNaN(iscomp)) {
            highlights.push({ seqIndex, startPosition, endPosition, iscomp });
        }
    });
    return highlights;
}

function generateHtml(sequences, highlights) {
    const uniqueIndices = [...new Set(highlights.map(h => h.seqIndex))];
    let htmlContent = '';

    uniqueIndices.forEach(index => {
        const sequence = sequences[index];
        const highlightsForSequence = highlights.filter(h => h.seqIndex === index);
        const highlightedSequence = highlightSequence(sequence.sequence, highlightsForSequence);

        htmlContent += `
            <div class="header">${sequence.header}</div>
            <div class="sequence">${highlightedSequence}</div>
        `;
    });

    return htmlContent;
}

function highlightSequence(sequence, highlights) {
    // Sort the highlights by startPosition
    highlights.sort((a, b) => a.startPosition - b.startPosition);

    // Merge overlapping intervals
    const mergedHighlights = [];
    highlights.forEach(({ startPosition, endPosition, iscomp }) => {
        if (
            mergedHighlights.length === 0 ||
            mergedHighlights[mergedHighlights.length - 1].endPosition < startPosition
        ) {
            mergedHighlights.push({ startPosition, endPosition, iscomp });
        } else {
            const last = mergedHighlights[mergedHighlights.length - 1];
            last.endPosition = Math.max(last.endPosition, endPosition);
            // If `iscomp` differs, we keep the last `iscomp` for simplicity
            if (last.iscomp !== iscomp) {
                last.iscomp = iscomp;
            }
        }
    });

    // Apply the highlights to the sequence
    let highlighted = '';
    let lastIndex = 0;
    
    mergedHighlights.forEach(({ startPosition, endPosition, iscomp }) => {
        if (startPosition > lastIndex) {
            highlighted += sequence.substring(lastIndex, startPosition);
        }

        // Apply a different color based on `iscomp`
        const colorClass = iscomp === 1 ? 'highlight-comp' : 'highlight';
        highlighted += `<span class="${colorClass}">${sequence.substring(startPosition, endPosition)}</span>`;
        lastIndex = endPosition;
    });

    if (lastIndex < sequence.length) {
        highlighted += sequence.substring(lastIndex);
    }


    return highlighted;
}
"""
