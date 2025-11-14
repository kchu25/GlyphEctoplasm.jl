# ============================================================================
# CSS TEMPLATE
# ============================================================================

"""
Complete CSS stylesheet for motif visualization pages.

Features:
- Responsive grid layout (auto-fit columns)
- Modal windows with overlay
- Interactive range sliders
- Navigation bar styling
- Syntax highlighting for sequences
- Hover effects and animations
"""
template_css = mt"""
body {
    font-family: Arial, sans-serif;
    background-color: #f5f5f7;
    margin: 0;
    padding: 0;
    justify-content: center;
    align-items: center;
    height: 100vh;
}

span.putBar {
  border-top: 1px solid #5144FA;;
}

.wrapper {
    width: 65%; /* Center 60% of the page */
    margin: 0 auto; /* Ensure it's centered */
    overflow: visible; /* Allow enlarged cards to extend beyond wrapper */
    padding: 0 60px; /* Add horizontal padding to prevent edge cards from being clipped */
}

.current {
    font-weight: bold;
    pointer-events: none;
    text-decoration: none; /* Remove underline */
}

.container {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
    gap: 15px;
    padding: 20px;
    justify-items: center;
    overflow: visible; /* Allow enlarged cards to extend beyond container */
    margin: 0 -60px; /* Extend container 60px beyond wrapper on both sides */
    padding: 20px 80px; /* Compensate with extra horizontal padding */
}

.sliderGroup {
    width: 100%;
    max-width: 150px;
    height: 60px;
    cursor: pointer;
    transition: transform 0.3s ease, box-shadow 0.3s ease, border-color 0.3s ease;
    border-radius: 6px;
    overflow: visible;
    background-color: #ffffff;
    /* Border and gradient shadow will be set dynamically via JavaScript */
    border: 1px solid #e5e5e5;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 2px 4px;
    position: relative;
}

.sliderGroup:hover {
    transform: scale(1.4);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    border-color: #8e8e93;
    border-width: 0.5px;
    z-index: 100;
}

.imageTextContainer {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 100%;
    height: 100%;
}

.imageContainer {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    position: relative;  /* Enable absolute positioning for overlay */
}

.imageContainer img {
    width: 100%;
    height: 100%;
    object-fit: contain;
    border-radius: 4px;
}

/* Filter index overlay for multi-motif cards */
.filter-index-overlay {
    position: absolute;
    top: 4px;
    right: 4px;
    background-color: rgba(211, 211, 211, 0.7);  /* Light gray with transparency */
    color: #666;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 10px;
    font-weight: 500;
    pointer-events: none;  /* Don't interfere with clicks */
    user-select: text;  /* Allow text selection for Ctrl+F */
    z-index: 10;
}

.textContainer {
    display: none;
}

/* Keep special styling for the current page link */
.textContainer a {
    text-decoration: none; /* Remove underline */
    color: #5144FA; /* Light blue color for links */
    margin: 0 5px; /* Even spacing between links */
}
    
.imageText {
    background-color: rgba(0, 0, 0, 0.05);
    padding: 6px;
    border-radius: 15px;
    font-size: 8px;
    white-space: nowrap; /* Ensure text doesn't wrap to fit the width */
    overflow: hidden; /* Hide any overflow */
    text-overflow: ellipsis; /* Add ellipsis if the text overflows */
}

.sliderContainer {
    display: none;
}

.sliderContainer input[type="range"] {
    width: 25%; /* Shortened width of the sliders */
    margin-top: 10px;
    -webkit-appearance: none;
    appearance: none;
    height: 10px;
    background: #ddd;
    outline: none;
    opacity: 0.65;
    /* transition: opacity 0.01s ease; */
    border: 1.5px solid gray; /* Thin black border around the slider */
    border-radius: 6px;
}

.sliderContainer input[type="range"]:hover {
    opacity: 1;
}

.sliderContainer input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    appearance: none;
    width: 17px;
    height: 17px;
    border-radius: 60%;    
    background: silver;
    border: 1.5px solid black;    
    cursor: pointer;
    /* transition: background 0.01s ease; */
    z-index: 1;
}

.sliderContainer input[type="range"]::-moz-range-thumb {
    width: 20px;
    height: 20px;
    border-radius: 50%;
    background: silver;
    border: 2px solid black;
    cursor: pointer;
    /* transition: background 0.01s ease; */
}

 /* Modal Styles (for sequence substring highlight)*/

.column {
    display: flex;
    flex-direction: column;
    width: 48%; /* Adjust width as needed */
}

.modal {
    display: none; /* Hidden by default */
    position: fixed; /* Stay in place */
    z-index: 1000; /* Sit on top */
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    overflow: auto; /* Enable scroll if needed */
    background-color: rgba(0,0,0,0.4); /* Black with opacity */
}

.modal-content {
    background-color: #fefefe;
    margin: 20% auto; /* 20% from the top and centered */
    padding: 25px;
    border: 2px solid #888;
    width: 80%; /* Could be more or less, depending on screen size */
    max-width: 985px; /* Limit max width */
    font-family: Arial, sans-serif;
    position: relative;
}

.close {
    color: #aaa;
    float: right;
    font-size: 28px;
    font-weight: bold;
    cursor: pointer;
}

.close:hover,
.close:focus {
    color: black;
    text-decoration: none;
    cursor: pointer;
}

.highlight {
    font-weight: bold;
    color: orange;
}

.highlight-comp {
    font-weight: bold;
    color: LightSteelBlue;
}

.sequence {
    font-family: monospace;
    white-space: pre-wrap; /* Preserve whitespace and wrap text */
    margin: 0;
    padding: 5px;
    text-align: left; /* Ensure text is aligned to the left */
}

.header {
    font-family: monospace;
    margin: 1px 0; /* Space above and below headers */
}

/* Modal cluster styles */
#highlightModal_cluster {
   display: none;
   position: fixed;
   z-index: 1;
   left: 0;
   top: 0;
   width: 100%;
   height: 100%;
   background-color: rgba(0, 0, 0, 0.5);
   overflow: auto;
}

#highlightModal_text {
   display: none;
   position: fixed;
   z-index: 1;
   left: 0;
   top: 0;
   width: 100%;
   height: 100%;
   background-color: rgba(0, 0, 0, 0.5);
   overflow: auto;

    /* Allow dynamic width adjustment */
   justify-content: center; 
   align-items: center; /* Center the modal */
}

#highlightModal_img {
   display: none;
   position: fixed;
   z-index: 1;
   left: 0;
   top: 0;
   width: 100%;
   height: 100%;
   background-color: rgba(0, 0, 0, 0.5);
   overflow: auto;
}


#highlightModal_text_content {
    background-color: white;
    margin: 15% auto;
    padding: 20px;
    border: 1px solid #888;
    text-align: center;
    max-width: 90%; /* Prevent it from exceeding 90% of the screen width */
    min-width: 200px; /* Ensure it doesn't shrink too much */
    word-wrap: break-word; /* Break long words if needed */
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2); /* Optional: add some style */
    border-radius: 8px; /* Optional: rounded corners */
    display: flex;
}


#highlightModal_img_content {
    background-color: white;
    margin: 15% auto;
    padding: 20px;
    border: 1px solid #888;
    text-align: center;
    max-width: 90%; /* Prevent it from exceeding 90% of the screen width */
    min-width: 200px; /* Ensure it doesn't shrink too much */
    word-wrap: break-word; /* Break long words if needed */
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2); /* Optional: add some style */
    border-radius: 8px; /* Optional: rounded corners */
    display: flex;
}


#highlightContent {
   background-color: white;
   margin: 15% auto;
   padding: 20px;
   border: 1px solid #888;
   width: 80%;
   text-align: center;
   display: flex;
   justify-content: space-around;
   align-items: center;
}

#modalText {
    font-size: 12px; /* Adjust this value as needed */
}

#modalText1 {
    font-size: 12px; /* Adjust this value as needed */
}

#copyButton {
    padding: 12px 24px;
    background-color: white;
    color: gray; /* gray text color */
    border: 2px solid lightgray; /* gray border */
    border-radius: 50px; /* Rounded, pill-like shape */
    cursor: pointer;
    font-size: 12px;
    font-weight: 500;
    margin-top: 15px;
    text-align: center;
    transition: all 0.3s ease; /* Smooth transition */
    box-shadow: 0 4px 10px rgba(0, 123, 255, 0.1); /* Subtle shadow */
}

#copyButton:hover {
    background-color: lightgray; /* Blue background on hover */
    color: black; /* White text on hover */
    border-color: lightgray; /* Darker border on hover */
    box-shadow: 0 6px 15px rgba(0, 123, 255, 0.2); /* Stronger shadow on hover */
}

#copyButton:focus {
    outline: none; /* Removes outline when focused */
    box-shadow: 0 0 5px rgba(0, 123, 255, 0.5); /* Focused glow effect */
}

.modal-column {
   flex: 1;
   padding: 10px;
}

.modal-column img {
   max-width: 100%;
}

.hover-window {
    position: fixed;
    left: 20px; /* Position from the left */
    top: 100px; /* Position from the top */
    width: 125px;
    padding: 10px;
    background-color: rgba(255, 255, 255, 0.9);
    border: 1px solid #ccc;
    border-radius: 5px;
    box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
    transition: transform 0.2s;
    font-size: 11px; /* Adjust this value as needed */
}

.hover-meta-data {    
    font-size: 8px; /* Adjust this value as needed */
}

.cl {
    list-style-type: none; /* Remove default list style */
    padding: 0; /* Remove padding */
    margin-bottom: 5px;
    display: flex; /* Use flex to align items */
    align-items: center; /* Center vertically */
}

.color-square {
    width: 15px; /* Width of the color square */
    height: 15px; /* Height of the color square */
    margin-right: 10px; /* Space between square and text */
    border-radius: 3px; /* Optional: rounded corners */
}

/* General navigation styling */
#nav {
    font-family: 'Arial', sans-serif; /* Use a clean and widely available font */
    font-size: 14px; /* Readable font size */
    background-color: #f9f9f9; /* Light background for contrast */
    border: 1px solid #ddd; /* Subtle border */
    border-radius: 12px; /* Rounded corners for a modern look */
    text-align: center; /* Center align the text */
    padding: 5px 10px; /* Add some padding for spacing */
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1); /* Subtle shadow for depth */
    max-width: fit-content; /* Limit the width of the navigation box */
    margin: 0 auto; /* Center the nav box within the page */
}

/* Styling for navigation links */
#nav a {
    text-decoration: none; /* Remove underline */
    color: #007BFF; /* Light blue color for links */
    margin: 0 5px; /* Even spacing between links */
}

/* Current page link styling */
#nav a.current {
    font-weight: bold; /* Highlight the current page */
    color: #000; /* Darker color for the current page */
}

/* The horizontal line */
.horizontal-line {
   border-top: 1px solid #ddd; /* Thin horizontal line */
   margin: 10px 0; /* Space above and below the line */
 }

/* ============================================================================
   GRID TOGGLE BUTTON STYLES
   ============================================================================ */

.grid-toggle-container {
    margin: 20px 0;
    padding: 0;
    width: 100%;
}

.grid-toggle-bar {
    display: flex;
    align-items: center;
    background-color: transparent;
    padding: 8px 20px;
    transition: all 300ms ease;
}

.grid-toggle-button {
    background-color: #ffffff;
    border: none;
    padding: 10px 20px;
    font-size: 14px;
    font-weight: 500;
    color: #8e8e93;
    cursor: pointer;
    border-radius: 8px;
    transition: all 250ms ease;
    display: flex;
    align-items: center;
    gap: 8px;
    box-shadow: 0 2px 6px rgba(0, 0, 0, 0.08);
    border-bottom: 3px solid transparent;
}

.grid-toggle-button:hover {
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.12);
    transform: translateY(-1px);
}

.grid-toggle-button.active {
    color: #1d1d1f;
    border-bottom-color: #8e8e93;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}

.grid-toggle-icon {
    font-size: 14px;
    transition: transform 250ms ease;
}

.grid-toggle-button.active .grid-toggle-icon {
    transform: rotate(180deg);
}

.grid-wrapper {
    max-height: 0;
    overflow: hidden;
    transition: max-height 500ms cubic-bezier(0.4, 0, 0.2, 1), 
                transform 500ms cubic-bezier(0.4, 0, 0.2, 1),
                opacity 400ms ease;
    opacity: 0;
    transform: translateY(-20px);
}

.grid-wrapper.visible {
    max-height: 10000px;
    opacity: 1;
    transform: translateY(0);
    overflow: visible; /* Allow enlarged cards to extend beyond wrapper */
}

/* ============================================================================
   SINGLETON GRID STYLES
   ============================================================================ */

.singleton-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(120px, 1fr));
    gap: 15px;
    padding: 20px;
    justify-items: center;
}

.singleton-cell {
    width: 100%;
    max-width: 120px;
    height: 50px;
    cursor: pointer;
    transition: transform 0.3s ease, box-shadow 0.3s ease, border-color 0.3s ease;
    border-radius: 6px;
    overflow: visible;
    background-color: #ffffff;
    /* Border and gradient shadow will be set dynamically via JavaScript */
    border: 1px solid #e5e5e5;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 2px 4px;
    position: relative;
}

.singleton-cell:hover {
    transform: scale(1.4);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    border-color: #8e8e93;
    border-width: 0.5px;
    z-index: 100;
}

.singleton-img {
    width: 100%;
    height: 100%;
    object-fit: contain;
}

/* Filter index overlay for singleton cards */
.singleton-filter-overlay {
    position: absolute;
    top: 2px;
    right: 2px;
    background-color: rgba(211, 211, 211, 0.7);  /* Light gray with transparency */
    color: #666;
    padding: 1px 4px;
    border-radius: 2px;
    font-size: 9px;
    font-weight: 500;
    pointer-events: none;  /* Don't interfere with clicks */
    user-select: text;  /* Allow text selection for Ctrl+F */
    z-index: 10;
}

/* ============================================================================
   SINGLETON MODAL STYLES
   ============================================================================ */

.singleton-modal {
    display: none;
    position: fixed;
    z-index: 2000;
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    overflow: auto;
    background-color: rgba(0, 0, 0, 0.6);
    animation: fadeIn 0.3s;
}

@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}

.singleton-modal-content {
    background-color: #f9f9f9;
    margin: 5% auto;
    padding: 0;
    border-radius: 12px;
    width: 90%;
    max-width: 1100px;
    box-shadow: 0 10px 40px rgba(0, 0, 0, 0.3);
    animation: slideDown 0.3s;
    position: relative;
}

@keyframes slideDown {
    from {
        transform: translateY(-50px);
        opacity: 0;
    }
    to {
        transform: translateY(0);
        opacity: 1;
    }
}

.singleton-close {
    color: #aaa;
    position: absolute;
    top: 15px;
    right: 20px;
    font-size: 32px;
    font-weight: bold;
    cursor: pointer;
    z-index: 2001;
    transition: color 0.2s;
}

.singleton-close:hover,
.singleton-close:focus {
    color: #000;
}

.singleton-modal-body {
    padding: 40px 30px 30px 30px;
    display: flex;
    flex-direction: row;
    align-items: flex-start;
    gap: 25px;
}

.singleton-modal-left {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 20px;
    padding: 20px;
}

.singleton-modal-influence-container {
    width: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 10px;
    background-color: #ffffff;
    padding: 1rem;
    border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.singleton-modal-influence-container img {
    max-width: 100%;
    max-height: 500px;
    object-fit: contain;
}

.singleton-modal-influence-label {
    font-size: 13px;
    color: #555;
    font-weight: 500;
    text-align: center;
}

.singleton-modal-right {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    gap: 20px;
}

.singleton-modal-right h3 {
    margin: 0;
    color: #333;
    font-size: 20px;
    text-align: center;
    padding-top: 10px;
}

.singleton-modal-img-container {
    width: 100%;
    max-width: 500px;
    height: 180px;
    padding: 1rem;
    border-radius: 8px;
    display: flex;
    justify-content: center;
    align-items: center;
    background-color: #ffffff;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.singleton-modal-img-container img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
    border-radius: 4px;
}

.singleton-modal-info {
    width: 100%;
    padding: 0 10px;
}

.singleton-info-item {
    background-color: #ffffff;
    padding: 12px 15px;
    margin: 8px 0;
    border-radius: 8px;
    font-size: 14px;
    color: #444;
    border-left: 3px solid #8e8e93;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.singleton-info-item a {
    color: #007aff;
    text-decoration: none;
    font-weight: 500;
}

.singleton-info-item a:hover {
    text-decoration: underline;
}

.singleton-info-item:empty {
    display: none;
}

/* ============================================================================
   MULTI-MOTIF MODAL STYLES
   ============================================================================ */

.multi-modal {
    display: none;
    position: fixed;
    z-index: 2000;
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    overflow: auto;
    background-color: rgba(0, 0, 0, 0.6);
    animation: fadeIn 0.3s;
}

@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}

.multi-modal-content {
    background-color: #f9f9f9;
    margin: 5% auto;
    padding: 0;
    border-radius: 12px;
    width: 90%;
    max-width: 1100px;
    box-shadow: 0 10px 40px rgba(0, 0, 0, 0.3);
    animation: slideDown 0.3s;
    position: relative;
}

@keyframes slideDown {
    from {
        transform: translateY(-50px);
        opacity: 0;
    }
    to {
        transform: translateY(0);
        opacity: 1;
    }
}

.multi-close {
    color: #aaa;
    position: absolute;
    top: 15px;
    right: 20px;
    font-size: 32px;
    font-weight: bold;
    cursor: pointer;
    z-index: 2001;
    transition: color 0.2s;
}

.multi-close:hover,
.multi-close:focus {
    color: #000;
}

.multi-modal-body {
    padding: 40px 30px 30px 30px;
    display: flex;
    flex-direction: row;
    align-items: flex-start;
    gap: 40px;
}

.multi-modal-left {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 25px;
}

.multi-modal-right {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-start;
    gap: 30px;
}

.multi-modal-right h3 {
    margin: 0;
    color: #333;
    font-size: 20px;
    text-align: center;
    padding-top: 10px;
}

.multi-modal-influence-container {
    width: 100%;
    padding: 1rem;
    border-radius: 8px;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 10px;
    background-color: #ffffff;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.multi-modal-influence-container img {
    max-width: 100%;
    max-height: 200px;
    object-fit: contain;
}

.multi-modal-influence-label {
    font-size: 13px;
    color: #555;
    font-weight: 500;
    text-align: center;
}

.multi-modal-img-container {
    width: 100%;
    max-width: 450px;
    height: 140px;
    padding: 1rem;
    border-radius: 8px;
    display: flex;
    justify-content: center;
    align-items: center;
    background-color: #ffffff;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.multi-modal-img-container img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
    border-radius: 4px;
}

.multi-modal-slider {
    width: 100%;
    max-width: 600px;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 10px;
    padding: 1rem;
    border-radius: 8px;
    background-color: #ffffff;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.multi-modal-slider input[type="range"] {
    width: 80%;
    -webkit-appearance: none;
    appearance: none;
    height: 8px;
    background: #ddd;
    outline: none;
    border-radius: 4px;
    border: 1px solid #ccc;
}

.multi-modal-slider input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    appearance: none;
    width: 18px;
    height: 18px;
    border-radius: 50%;
    background: #8e8e93;
    border: 2px solid #ffffff;
    cursor: pointer;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

.multi-modal-slider input[type="range"]::-moz-range-thumb {
    width: 18px;
    height: 18px;
    border-radius: 50%;
    background: #8e8e93;
    border: 2px solid #ffffff;
    cursor: pointer;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
}

.multi-modal-slider-label {
    font-size: 18px;
    color: #333;
    font-weight: 700;
}

.multi-modal-info {
    width: 100%;
    padding: 0 10px;
}

.multi-info-item {
    background-color: #ffffff;
    padding: 12px 15px;
    margin: 8px 0;
    border-radius: 8px;
    font-size: 14px;
    color: #444;
    border-left: 3px solid #8e8e93;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.multi-info-item a {
    color: #007aff;
    text-decoration: none;
    font-weight: 500;
}

.multi-info-item a:hover {
    text-decoration: underline;
}

.multi-info-item:empty {
    display: none;
}


"""
