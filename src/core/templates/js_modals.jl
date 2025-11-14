# ============================================================================
# JAVASCRIPT MODAL FUNCTIONS
# ============================================================================

"""
JavaScript functions for managing various modal windows.

Modal types:
1. Image modal - Display images in fullscreen
2. Text modal - Display copyable text content
3. Cluster modal - Display image + text side-by-side
4. Singleton modal - Display singleton motif details
5. Multi-motif modal - Display multi-motif details with sliders

Functions:
- `openHtmlWindowImg()`: Open image modal
- `openHtmlWindowText()`: Open text modal with copy functionality
- `openHtmlWindow()`: Open cluster modal (image + text)
- `openSingletonModal()`: Open singleton motif modal
- `openMultiModal()`: Open multi-motif modal
- `closeModal*()`: Close specific modals
- `copyText()`: Copy text to clipboard
- Modal event handlers for clicks and keyboard shortcuts
"""
script_modals_str = raw"""
function openHtmlWindowImg(imageFile) {
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;
    // Set the image source in the modal
    const modalImage = document.getElementById('modalImage1');
    modalImage.src = imageFile;

    // Wait for the image to load to get its natural dimensions
    modalImage.onload = function () {
        // Get the image's natural width
        const imgWidth = modalImage.naturalWidth;

        // Dynamically set the modal width (optional max width for responsiveness)
        const modal = document.getElementById('highlightModal_img_content');
        modal.style.width = imgWidth > 800 ? '800px' : imgWidth + 'px'; // Cap at 800px for responsiveness
    };

    // Display the modal
    document.getElementById('highlightModal_img').style.display = "block";
}
   
function openHtmlWindowText(textContent) {
    // Store the current scroll position
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;

    // Set the text content
    const modalText = document.getElementById('modalText1');
    modalText.innerHTML = textContent;

    // Dynamically adjust the modal width based on text length
    const modalContent = document.getElementById('highlightModal_text_content');
    const contentLength = textContent.length;

    // Calculate width: 10px per character, with min and max limits
    const width = Math.min(Math.max(contentLength * 10, 200), 800); 
    modalContent.style.width = width + 'px';

    // Display the modal
    document.getElementById('highlightModal_text').style.display = "flex";
}

// Function to copy text from the modal
function copyText() {
    const modalTextElement = document.getElementById('modalText1');
    const originalText = modalTextElement.innerText;  // Save the original text
    const textToCopy = originalText;
    const textarea = document.createElement('textarea');
    textarea.value = textToCopy;
    document.body.appendChild(textarea);
    textarea.select();
    document.execCommand('copy');
    document.body.removeChild(textarea);

    // Display success message temporarily in the modal
    modalTextElement.innerHTML = 'string copied successfully!';

    // After 1 second, revert the content back to the original
    setTimeout(() => {
        modalTextElement.innerHTML = originalText;
    }, 1000);  // 1 second delay
}

function openHtmlWindow(imageFile, textContent) {
    // Store the current scroll position
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;

    // Set the image source in the modal
    document.getElementById('modalImage').src = imageFile;

    // Set the dynamic text content
    document.getElementById('modalText').innerHTML = textContent;

    // Display the modal
    document.getElementById('highlightModal_cluster').style.display = "block";
}

function closeModal() {
    document.getElementById('highlightModal').style.display = "none";

    // Restore the scroll position
    window.scrollTo(0, scrollPosition);
}

function closeModal_text() {
    // Hide the modal
    document.getElementById('highlightModal_text').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}

function closeModal_cluster() {
    // Hide the modal
    document.getElementById('highlightModal_cluster').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}

function closeModal_img() {
    // Hide the modal
    document.getElementById('highlightModal_img').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}

function closeMultiMotifModal() {
    document.getElementById('multiMotifModal').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}

// ============================================================================
// SINGLETON MODAL FUNCTIONS
// ============================================================================

let singletonScrollPosition = 0;

function openSingletonModal(index) {
    // Store scroll position
    singletonScrollPosition = window.pageYOffset || document.documentElement.scrollTop;
    
    // Get data from global array (1-indexed from Julia)
    const data = window.singletonData[index - 1];
    if (!data) return;
    
    // Populate modal
    document.getElementById('singletonModalImg').src = data.img;
    document.getElementById('singletonModalImg').alt = data.alt;
    document.getElementById('singletonModalInfluence').src = data.influence;
    document.getElementById('singletonModalTitle').innerHTML = data.alt;
    document.getElementById('singletonModalText1').innerHTML = data.text1 || '';
    document.getElementById('singletonModalText2').innerHTML = data.text2 || '';
    document.getElementById('singletonModalText3').innerHTML = data.text3 || '';
    document.getElementById('singletonModalText4').innerHTML = data.text4 || '';
    document.getElementById('singletonModalText5').innerHTML = data.text5 || '';
    
    // Show modal
    document.getElementById('singletonModal').style.display = 'block';
}

function closeSingletonModal() {
    document.getElementById('singletonModal').style.display = 'none';
    // Restore scroll position
    window.scrollTo(0, singletonScrollPosition);
}

// ============================================================================
// MULTI-MOTIF MODAL FUNCTIONS
// ============================================================================

let multiModalScrollPosition = 0;
let multiModalData = {
    images: [],
    labels: [],
    texts: [],
    baseFolders: []  // Store base folder paths for influence images
};

function openMultiModal(index) {
    // Store scroll position
    multiModalScrollPosition = window.pageYOffset || document.documentElement.scrollTop;
    
    // Get data for this motif (1-indexed from display, map to array index)
    let modeIndex = index - 1;  // Convert from 1-indexed to 0-indexed
    
    // Check if we have an index mapping (for unified template)
    if (typeof window.multiMotifIndexMap !== 'undefined' && window.multiMotifIndexMap[index] !== undefined) {
        modeIndex = window.multiMotifIndexMap[index];
    }
    
    if (modeIndex < 0 || modeIndex >= multiModalData.images.length) {
        console.error('Invalid multi-motif index:', index, 'mapped to:', modeIndex);
        return;
    }
    
    const images = multiModalData.images[modeIndex];
    const labels = multiModalData.labels[modeIndex];
    const texts = multiModalData.texts[modeIndex];
    const baseFolder = multiModalData.baseFolders[modeIndex];
    
    // Initialize modal with first image
    updateMultiModalContent(0, images, labels, texts, baseFolder);
    
    // Setup slider
    const slider = document.getElementById('multiMotifSlider');
    slider.max = images.length - 1;
    slider.value = 0;
    
    // Hide slider if only one image
    const sliderContainer = document.querySelector('.multi-modal-slider');
    if (images.length === 1) {
        sliderContainer.style.display = 'none';
    } else {
        sliderContainer.style.display = 'flex';
    }
    
    // Attach slider event
    slider.oninput = function() {
        updateMultiModalContent(parseInt(this.value), images, labels, texts, baseFolder);
    };
    
    // Show modal
    document.getElementById('multiMotifModal').style.display = 'block';
}

function updateMultiModalContent(sliderValue, images, labels, texts, baseFolder) {
    // Update logo image
    document.getElementById('multiMotifImage').src = images[sliderValue];
    
    // Format title for slider label: convert "pattern 8 and 112 are 0 nucleotides apart<br>" 
    // to "pattern 8 and 112: 0 nucleotides apart" (keep <br> for line breaks in triplets+)
    let titleText = labels[sliderValue].trim();
    titleText = titleText.replace(/ are /g, ': ');
    
    document.getElementById('multiMotifText1').innerHTML = texts[sliderValue][0] || '';
    document.getElementById('multiMotifText2').innerHTML = texts[sliderValue][1] || '';
    document.getElementById('multiMotifText3').innerHTML = texts[sliderValue][2] || '';
    document.getElementById('multiMotifText4').innerHTML = texts[sliderValue][3] || '';
    document.getElementById('multiMotifText5').innerHTML = texts[sliderValue][4] || '';
    document.getElementById('multiMotifRangeLabel').innerHTML = titleText;
    
    // Update fixed distance influence plot (changes with slider)
    const fixedInfluencePath = images[sliderValue].replace('.png', '_influence.png');
    document.getElementById('multiMotifInfluenceFixed').src = fixedInfluencePath;
    
    // Update relaxed distance influence plot (constant for all sliders in this motif)
    const relaxedInfluencePath = baseFolder + '/influence_relaxed.png';
    document.getElementById('multiMotifInfluenceRelaxed').src = relaxedInfluencePath;
}

function closeMultiModal() {
    document.getElementById('multiMotifModal').style.display = 'none';
    // Restore scroll position
    window.scrollTo(0, multiModalScrollPosition);
}

// Alias for compatibility with HTML onclick
function closeMultiMotifModal() {
    closeMultiModal();
}

// Store multi-motif data when loaded
function storeMultiModalData(images, labels, texts) {
    multiModalData.images = images;
    multiModalData.labels = labels;
    multiModalData.texts = texts;
    
    // Extract base folder paths from first image of each motif
    multiModalData.baseFolders = images.map(modeImages => {
        if (modeImages.length > 0) {
            // Get folder path from first image: "pair_motifs/2_112/0.png" -> "pair_motifs/2_112"
            const firstImage = modeImages[0];
            const lastSlashIndex = firstImage.lastIndexOf('/');
            return lastSlashIndex >= 0 ? firstImage.substring(0, lastSlashIndex) : '';
        }
        return '';
    });
}

// Close the modal when clicking outside of it
window.onclick = function(event) {
    const modal = document.getElementById('highlightModal');
    if (event.target === modal) {
        closeModal();
    }        
    const modal_cluster = document.getElementById('highlightModal_cluster');
    if (event.target === modal_cluster) {
        closeModal_cluster();
    }
    const modal_text = document.getElementById('highlightModal_text');
    if (event.target === modal_text) {
        closeModal_text();
    }
    const modal_img = document.getElementById('highlightModal_img');
    if (event.target === modal_img) {
        closeModal_img();
    }
    const multiModal = document.getElementById('multiMotifModal');
    if (multiModal && event.target === multiModal) {
        closeMultiMotifModal();
    }
    const singletonModal = document.getElementById('singletonModal');
    if (singletonModal && event.target === singletonModal) {
        closeSingletonModal();
    }
}

// Close the modal when pressing the "Esc" key
window.onkeydown = function(event) {
    if (event.key === "Escape") {
        closeModal();
        closeModal_cluster();
        closeModal_text();
        closeModal_img();
        closeMultiMotifModal();
        const singletonModal = document.getElementById('singletonModal');
        if (singletonModal) closeSingletonModal();
    }
}
"""
