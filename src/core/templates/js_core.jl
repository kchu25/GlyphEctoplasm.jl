# ============================================================================
# JAVASCRIPT CORE FUNCTIONS
# ============================================================================

"""
Core JavaScript functions for navigation, data loading, and slider initialization.

Mustache variables:
- `j`: Current page number
- `upto`: Total number of pages
- `mode_counts`: Number of motif modes to display

Functions:
- `createArray()`: Helper to generate page number arrays
- `updateNav()`: Generate navigation links
- `loadData()`: Load motif data from JSON
- `initializeSliders()`: Set up range sliders for multi-motif cards
"""
script_core_str = raw"""
// Function to load images, labels, and texts from JSON

function createArray(num) {
    return Array.from({ length: num }, (_, i) => i + 1);
}

function updateNav() {
    const currentPage = 'index{{:j}}.html'; // Set this to the current page name
    const availablePages = createArray({{:upto}}); // List of available page numbers

    const navHtml = availablePages.map(num => {
        const page = `index${num}.html`;
        const label = num === 1 ? 'Pattern influence' :
                      num === 2 ? 'Generalization' :
                      num === 3 ? 'Statistics' :
                      num === 4 ? 'Readme' :
                      num === 5 ? 'Page 5' :
                      'Page 6';

        return `<a href="${page}" ${currentPage === page ? 'class="current"' : ''}>${label}</a>`;
    }).join(' &nbsp&nbsp | &nbsp&nbsp ');

    // Note: Readme link removed since page 3 is now "Readme"
    document.getElementById('nav').innerHTML = navHtml;
}

window.onload = updateNav;

function loadData(modeCount, jsonData, callback) {
    const images = [];
    const labels = [];
    const texts = [];

    // Try to detect if we're using grouped mode names (e.g., mode_g1_1) or simple ones (mode_1)
    const firstModeKey = Object.keys(jsonData).find(k => k.startsWith('mode_'));
    
    if (firstModeKey) {
        // Get all mode keys and sort them
        const modeKeys = Object.keys(jsonData)
            .filter(k => k.startsWith('mode_'))
            .sort((a, b) => {
                // Extract the final number from mode strings
                const aMatch = a.match(/_(\d+)$/);
                const bMatch = b.match(/_(\d+)$/);
                if (aMatch && bMatch) {
                    return parseInt(aMatch[1]) - parseInt(bMatch[1]);
                }
                return a.localeCompare(b);
            });
        
        // Use sorted keys to build arrays
        for (const modeKey of modeKeys) {
            const modeImages = jsonData[modeKey].pwms;
            const modeLabels = jsonData[modeKey].labels;
            const modeTexts = jsonData[modeKey].texts;

            images.push(modeImages);
            labels.push(modeLabels);
            texts.push(modeTexts);
        }
    } else {
        // Fallback: try numeric mode_1, mode_2, etc.
        for (let modeIndex = 1; modeIndex <= modeCount; modeIndex++) {
            const modeKey = `mode_${modeIndex}`;
            if (jsonData[modeKey]) {
                const modeImages = jsonData[modeKey].pwms;
                const modeLabels = jsonData[modeKey].labels;
                const modeTexts = jsonData[modeKey].texts;

                images.push(modeImages);
                labels.push(modeLabels);
                texts.push(modeTexts);
            }
        }
    }

    callback(images, labels, texts);
}

// Initialize the sliders with the loaded data
function initializeSliders(images, labels, texts) {
    images.forEach((modeImages, index) => {
        const modeLabels = labels[index];
        const modeTexts = texts[index];
        const sliderId = `valR${index + 1}`;
        const imgElementId = `img${index + 1}`;
        const rangeElementId = `range${index + 1}`;
        const textElementId1 = `text${index + 1}_1`;
        const textElementId2 = `text${index + 1}_2`;
        const textElementId3 = `text${index + 1}_3`;
        const textElementId4 = `text${index + 1}_4`;
        const textElementId5 = `text${index + 1}_5`;
        const textElementId6 = `text${index + 1}_6`;
        const sliderElement = document.getElementById(sliderId);
        
        function showVal(newVal) {
            const imgElement = document.getElementById(imgElementId);
            imgElement.style.opacity = 0;  // Start fading out

            setTimeout(function() {
                document.getElementById(rangeElementId).innerHTML = modeLabels[newVal];
                document.getElementById(textElementId1).innerHTML = modeTexts[newVal][0];
                document.getElementById(textElementId2).innerHTML = modeTexts[newVal][1];
                document.getElementById(textElementId3).innerHTML = modeTexts[newVal][2];
                document.getElementById(textElementId4).innerHTML = modeTexts[newVal][3];
                document.getElementById(textElementId5).innerHTML = modeTexts[newVal][4];
                document.getElementById(textElementId6).innerHTML = modeTexts[newVal][5];
                imgElement.src = modeImages[newVal];
                imgElement.style.opacity = 1;  // Fade back in
            }, 250);  // Wait for half the transition duration (0.5s / 2)
        }

        // Set initial values for each slider
        document.getElementById(rangeElementId).innerHTML = modeLabels[0];
        document.getElementById(textElementId1).innerHTML = modeTexts[0][0];
        document.getElementById(textElementId2).innerHTML = modeTexts[0][1];
        document.getElementById(textElementId3).innerHTML = modeTexts[0][2];
        document.getElementById(textElementId4).innerHTML = modeTexts[0][3];
        document.getElementById(textElementId5).innerHTML = modeTexts[0][4];
        document.getElementById(textElementId6).innerHTML = modeTexts[0][5];
        document.getElementById(imgElementId).src = modeImages[0];

        // Attach event listener to each slider
        document.getElementById(sliderId).addEventListener('input', function() {
            showVal(this.value);
        });

        if (modeImages.length == 1) {
            sliderElement.style.display = 'none';
            document.getElementById(rangeElementId).innerHTML = modeLabels[0];
        }

    });
}

// Apply dynamic styling when the page loads and initialize data
document.addEventListener('DOMContentLoaded', function() {
    // Load data and initialize sliders
    loadData({{:mode_counts}}, jsonData, (images, labels, texts) => {
        initializeSliders(images, labels, texts);
        storeMultiModalData(images, labels, texts);
    });
    
    // Apply dynamic card styling
    applyDynamicCardStyling();
});

// hover window

    const hoverWindow = document.getElementById('hoverWindow');

    document.addEventListener('scroll', () => {
        const scrollY = window.scrollY || window.pageYOffset;
    });
"""
