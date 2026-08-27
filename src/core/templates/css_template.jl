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
/* ── Low-generalization warning ──────────────────────────────────────────
   Shown on the motif and summary pages when held-out R2 is below threshold.
   The motif list is a ranking (top/bottom N by median Banzhaf), not a
   significance test, so it reports groups even on data with no signal. */
.transform-note{
    display:flex; gap:0.75rem; align-items:flex-start;
    margin:0 auto 1.1rem auto; max-width:1100px;
    padding:0.7rem 0.95rem; border-radius:8px;
    background:#eef4fa; color:#173a5e; border:1px solid #9fc0dd;
    font-size:0.92rem; line-height:1.45;
}
.transform-note-badge{
    flex:0 0 auto; font-weight:700; letter-spacing:0.02em;
    text-transform:uppercase; font-size:0.74rem;
    padding:0.16rem 0.5rem; border-radius:5px;
    background:#9fc0dd; color:#123;
}
.transform-note-body{flex:1 1 auto;}
/* Footnote placement: separated from the content above, and quieter, since by
   the time a reader reaches it they have already seen the numbers it describes. */
.transform-note-bottom{
    margin:2.25rem auto 0.5rem auto;
    font-size:0.86rem;
    border-top-width:1px;
}
@media (prefers-color-scheme: dark){
    .transform-note{background:#152230; color:#cfe2f4; border-color:#3f6f9c;}
    .transform-note-badge{background:#3f6f9c; color:#eaf3fb;}
}
.generalization-warning{
    display:flex; gap:14px; align-items:flex-start;
    max-width:1100px; margin:18px auto 6px auto; padding:14px 18px;
    border:1px solid #d98324; border-left:6px solid #d98324; border-radius:8px;
    background:#fff8ef; color:#5a3a12;
    font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
    font-size:14px; line-height:1.5;
}
.generalization-warning-badge{
    flex:0 0 auto; align-self:center;
    background:#d98324; color:#fff; font-weight:700; font-size:11px;
    letter-spacing:.04em; text-transform:uppercase;
    padding:5px 10px; border-radius:4px; white-space:nowrap;
}
.generalization-warning-body{flex:1 1 auto;}
.generalization-warning-body strong{color:#8a4b00;}
.generalization-warning-body code{
    background:#f2e2cd; padding:1px 5px; border-radius:3px;
    font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:13px;
}
@media (prefers-color-scheme: dark){
    .generalization-warning{background:#2b2115; color:#f0dcc0; border-color:#c9791f;}
    .generalization-warning-body strong{color:#ffc478;}
    .generalization-warning-body code{background:#3d2f1d;}
}

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

/* Interaction indicator for pair-motifs with detected interactions */
.interaction-indicator {
    position: absolute;
    bottom: 1px;
    right: 6px;
    font-size: 21px;
    font-weight: 700;
    font-family: monospace;
    color: rgba(128, 128, 128, 0.4);  /* More transparent */
    pointer-events: none;
    user-select: none;
    z-index: 5;
    letter-spacing: -2px;
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


/* Insignificant indicator for singleton cards that failed significance filter */
.insignificant-indicator {
    position: absolute;
    bottom: -0.85px;
    right: 1.75px;
    font-size: 7.5px;
    font-weight: 400;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: #4D4847;
    letter-spacing: 0.2px;
    text-transform: lowercase;
    pointer-events: none;
    user-select: none;
    z-index: 5;
}


/* ============================================================================
   TOP-MOVERS SUMMARY PAGE STYLES
   ============================================================================ */

/* Summary page uses more of the viewport than the grouped page's wrapper. */
.top-movers-wrapper {
    width: 90%;
    max-width: 1130px;          /* match the narrower section so the right void shrinks */
    padding: 0 24px;
}

/* Slick protein masthead at the top of the summary page — single line. */
.protein-title {
    display: flex;
    flex-direction: row;
    flex-wrap: wrap;
    align-items: baseline;
    justify-content: center;
    gap: 8px;
    margin: 18px auto 6px auto;
}

.protein-title-eyebrow {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 10px;
    font-weight: 600;
    letter-spacing: 2px;
    text-transform: uppercase;
    color: #9ca3af;
}

.protein-title-name {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 22px;
    font-weight: 300;
    letter-spacing: 0.5px;
    color: #1f2937;
}

.protein-title-sep {
    color: #d1d5db;
    font-size: 14px;
}

.protein-title-meta {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 11px;
    color: #9ca3af;
    font-variant-numeric: tabular-nums;
}

.top-movers-section {
    max-width: 1080px;          /* fit the content; trims empty space on the right */
    margin: 0 auto 32px auto;
}

/* First section sits close under the protein masthead (no big <br><br> gap). */
.top-movers-section:first-of-type {
    margin-top: 8px;
}

.top-movers-heading {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 16px;
    font-weight: 600;
    letter-spacing: 0.2px;
    padding-bottom: 6px;
    border-bottom: 1px solid #e5e5e5;
}

.top-movers-heading-pos { color: #8b0000; }
.top-movers-heading-neg { color: #00008b; }

.top-movers-list {
    display: flex;
    flex-direction: column;
    gap: 14px;
    padding: 14px 0;
}

/* Fixed grid tracks so the three columns line up vertically across every row.
   Track 2 is sized to fit the name/epistasis so the WT track starts right after
   it (no big gap); track 3 takes the remaining width. */
.top-mover-row {
    display: grid;
    grid-template-columns: 190px 200px 1fr;   /* card track (190) > card image (150) widens card→epistasis gap */
    align-items: center;
    column-gap: 28px;
}

/* Dual-view rows (mutagenesis): two logo cards — reduced view + region view —
   ahead of the meta and wild-type tracks. */
.top-mover-row.dual {
    grid-template-columns: 190px 190px 200px 1fr;
}

/* Column header row above a dual-view list. */
.top-mover-header {
    align-items: end;
    margin-bottom: 2px;
}

.top-mover-colhead {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    color: #6b7280;
    text-align: center;
    width: 150px;                    /* match the card width so labels center over cards */
}

/* Left: clickable card. Reuses the singleton card look; border is set via JS. */
.top-mover-card {
    width: 150px;
    height: 60px;
    cursor: pointer;
    transition: transform 0.2s ease, box-shadow 0.2s ease;
    border-radius: 6px;
    overflow: visible;
    background-color: #ffffff;
    border: 1px solid #e5e5e5;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 2px 4px;
    position: relative;
}

.top-mover-card:hover {
    transform: scale(1.08);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.12);
    z-index: 50;
}

/* Middle: epistasis column (grid track 2). */
.top-mover-meta {
    min-width: 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: #1f2937;
}

/* Right: wild-type track (grid track 3). */
.top-mover-wt {
    min-width: 0;
}

.top-mover-name {
    font-size: 13px;
    font-weight: 600;
    margin-bottom: 4px;
}

/* Muted "position(s)" prefix before a mutation-span name. */
.top-mover-name-label {
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 0.03em;
    text-transform: uppercase;
    color: #9ca3af;
    margin-right: 5px;
}

/* One-click wild-type copy button (sits beside the protein title). */
.wt-copy-btn {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 11px;
    font-weight: 500;
    color: #475569;
    background: #ffffff;
    border: 1px solid #d8dee9;
    border-radius: 6px;
    padding: 4px 10px;
    cursor: pointer;
    transition: background 0.15s ease, border-color 0.15s ease, color 0.15s ease;
}

.wt-copy-btn:hover {
    background: #f1f5f9;
    border-color: #94a3b8;
}

.wt-copy-btn.copied {
    color: #166534;
    background: #f0fdf4;
    border-color: #86efac;
}

/* Epistasis (interaction coefficient) — two compact dark-gray lines per row. */
.top-mover-epistasis {
    display: flex;
    flex-direction: column;
    gap: 1px;
    margin-top: 2px;
}

.epi-line {
    font-size: 12px;
    color: #4b5563;                 /* dark gray */
    font-variant-numeric: tabular-nums;
}

.epi-line strong {
    color: #374151;                 /* slightly darker for the coefficient */
    font-weight: 600;
}

.epi-sub {
    font-size: 11px;
    color: #6b7280;
}

.epi-none { color: #9ca3af; }

.top-mover-empty {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: #9ca3af;
    font-size: 13px;
    padding: 8px 0;
}

/* Wild-type amino-acid track: residues with up-arrows under mutated positions. */
.wt-block { margin: 5px 0 9px 0; }
.wt-region { margin-bottom: 6px; }

.wt-span {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 10px;
    letter-spacing: 0.3px;
    color: #9ca3af;
    margin-bottom: 1px;
}

.wt-track {
    display: flex;
    flex-wrap: wrap;
    gap: 1px;
}

.wt-col {
    display: flex;
    flex-direction: column;
    align-items: center;
    min-width: 12px;
}

.wt-aa {
    font-family: 'SF Mono', ui-monospace, Menlo, Consolas, monospace;
    font-size: 11px;
    line-height: 1.2;
    color: #b6bcc6;          /* context residues recede */
}

.wt-col.is-mut .wt-aa {
    color: #111827;          /* mutated residue stands out */
    font-weight: 700;
}

.wt-arrow {
    height: 9px;
    line-height: 1;
    font-size: 8px;
    color: #d97706;          /* sign-neutral amber accent */
}

.wt-pos {
    line-height: 1;
    font-size: 7px;
    color: #d97706;
    font-variant-numeric: tabular-nums;
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
    background-color: rgba(15, 23, 42, 0.32);
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
    animation: fadeIn 0.3s;
}

@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}

.singleton-modal-content {
    background-color: #f5f7fa;
    margin: 5% auto;
    padding: 0;
    border-radius: 16px;
    width: 90%;
    max-width: 1100px;
    border: 1px solid rgba(15, 23, 42, 0.08);
    box-shadow:
        0 30px 80px -24px rgba(15, 23, 42, 0.28),
        0 8px 24px -12px rgba(15, 23, 42, 0.14);
    animation: slideDown 0.3s;
    position: relative;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'SF Pro Text', Roboto, sans-serif;
    -webkit-font-smoothing: antialiased;
    font-feature-settings: "tnum" 1, "kern" 1;
    color: #1f2937;
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
    color: #94a3b8;
    position: absolute;
    top: 14px;
    right: 16px;
    width: 32px;
    height: 32px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    font-size: 22px;
    font-weight: 400;
    line-height: 1;
    cursor: pointer;
    z-index: 2001;
    border-radius: 50%;
    transition: color 0.18s ease, background-color 0.18s ease;
}

.singleton-close:hover,
.singleton-close:focus {
    color: #0f172a;
    background-color: rgba(15, 23, 42, 0.06);
}

.singleton-modal-body {
    padding: 44px 36px 32px 36px;
    display: flex;
    flex-direction: row;
    align-items: flex-start;
    gap: 28px;
}

.singleton-modal-left {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 20px;
    padding: 0;
}

.singleton-modal-influence-container {
    width: 100%;
    box-sizing: border-box;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 10px;
    background-color: #ffffff;
    padding: 14px 16px;
    border: 1px solid rgba(15, 23, 42, 0.07);
    border-radius: 12px;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}

.singleton-modal-influence-container img {
    max-width: 100%;
    max-height: 500px;
    object-fit: contain;
}

.singleton-modal-influence-label {
    font-size: 11px;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: #94a3b8;
    font-weight: 600;
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
    color: #0f172a;
    font-size: 15px;
    font-weight: 500;
    letter-spacing: 0.02em;
    text-align: center;
    padding-top: 0;
    order: 2;
}

.singleton-modal-kde-container {
    width: 100%;
    max-width: 360px;
    height: 320px;
    box-sizing: border-box;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    gap: 10px;
    background-color: #ffffff;
    padding: 14px 16px;
    border: 1px solid rgba(15, 23, 42, 0.07);
    border-radius: 12px;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
    order: 0;
}

.singleton-modal-kde-container img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
}

.singleton-modal-kde-label {
    font-size: 11px;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: #94a3b8;
    font-weight: 600;
    text-align: center;
    margin-bottom: 6px;
}

.singleton-modal-img-container {
    width: 100%;
    max-width: 360px;
    height: 105px;
    padding: 14px 16px;
    border-radius: 12px;
    display: flex;
    justify-content: center;
    align-items: center;
    background-color: #ffffff;
    border: 1px solid rgba(15, 23, 42, 0.07);
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
    order: 1;
}

.singleton-modal-img-container img {
    max-width: 80%;
    max-height: 80%;
    object-fit: contain;
    border-radius: 4px;
}

.singleton-modal-info {
    width: 100%;
    padding: 0 6px;
}

.singleton-info-item {
    background-color: #ffffff;
    box-sizing: border-box;
    padding: 10px 14px;
    margin: 6px 0;
    border-radius: 10px;
    font-size: 13px;
    line-height: 1.5;
    color: #334155;
    border: 1px solid rgba(15, 23, 42, 0.06);
    box-shadow: 0 1px 1px rgba(15, 23, 42, 0.03);
    transition: border-color 0.18s ease, box-shadow 0.18s ease;
}

.singleton-info-item:hover {
    border-color: rgba(15, 23, 42, 0.12);
    box-shadow: 0 2px 6px rgba(15, 23, 42, 0.06);
}

/* Plain-language interpretation block in the card popups. Deliberately set
   apart from the .singleton-info-item / .multi-info-item metadata rows: this is
   the one line that reads as a sentence rather than a statistic, so it gets an
   accent rule and a small caps label instead of another white card. Hidden by
   the JS (display:none) when a motif has no sentence. */
.motif-description {
    background-color: #f8fafc;
    box-sizing: border-box;
    padding: 12px 14px;
    margin: 10px 0 6px 0;
    border-radius: 10px;
    border: 1px solid rgba(15, 23, 42, 0.08);
    border-left: 3px solid #64748b;
}

.motif-description-label {
    display: block;
    font-size: 10px;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: #64748b;
    margin-bottom: 5px;
}

.motif-description-text {
    display: block;
    font-size: 13.5px;
    line-height: 1.55;
    color: #1e293b;
}

.singleton-info-item strong,
.singleton-info-item b {
    color: #0f172a;
    font-weight: 600;
}

.singleton-info-item a {
    color: #2563eb;
    text-decoration: none;
    font-weight: 500;
    padding: 2px 8px;
    border-radius: 999px;
    background-color: rgba(37, 99, 235, 0.08);
    transition: background-color 0.18s ease, color 0.18s ease;
}

.singleton-info-item a:hover {
    background-color: rgba(37, 99, 235, 0.16);
    text-decoration: none;
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
    background-color: rgba(15, 23, 42, 0.32);
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
    animation: fadeIn 0.3s;
}

@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}

.multi-modal-content {
    background-color: #f5f7fa;
    margin: 5% auto;
    padding: 0;
    border-radius: 16px;
    width: 90%;
    max-width: 1100px;
    border: 1px solid rgba(15, 23, 42, 0.08);
    box-shadow:
        0 30px 80px -24px rgba(15, 23, 42, 0.28),
        0 8px 24px -12px rgba(15, 23, 42, 0.14);
    animation: slideDown 0.3s;
    position: relative;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'SF Pro Text', Roboto, sans-serif;
    -webkit-font-smoothing: antialiased;
    font-feature-settings: "tnum" 1, "kern" 1;
    color: #1f2937;
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
    color: #94a3b8;
    position: absolute;
    top: 14px;
    right: 16px;
    width: 32px;
    height: 32px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    font-size: 22px;
    font-weight: 400;
    line-height: 1;
    cursor: pointer;
    z-index: 2001;
    border-radius: 50%;
    transition: color 0.18s ease, background-color 0.18s ease;
}

.multi-close:hover,
.multi-close:focus {
    color: #0f172a;
    background-color: rgba(15, 23, 42, 0.06);
}

.multi-modal-body {
    padding: 40px 30px 30px 30px;
    display: flex;
    flex-direction: row;
    align-items: flex-start;
    gap: 40px;
}

.multi-modal-left {
    flex: 0 0 480px;
    width: 480px;
    min-width: 480px;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 25px;
}

.multi-modal-right {
    flex: 0 0 480px;
    width: 480px;
    min-width: 480px;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-start;
    gap: 30px;
}

.multi-modal-right h3 {
    margin: 0;
    color: #0f172a;
    font-size: 15px;
    font-weight: 500;
    letter-spacing: 0.02em;
    text-align: center;
    padding-top: 10px;
}

.multi-modal-influence-container {
    width: 100%;
    box-sizing: border-box;
    padding: 14px 16px;
    border-radius: 12px;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 10px;
    background-color: #ffffff;
    border: 1px solid rgba(15, 23, 42, 0.07);
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}

.multi-modal-influence-container img {
    max-width: 100%;
    max-height: 200px;
    object-fit: contain;
}

.multi-modal-influence-label {
    font-size: 11px;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: #94a3b8;
    font-weight: 600;
    text-align: center;
}

.multi-modal-kde-container {
    max-width: 480px;
    width: 480px;
    min-width: 480px;
    height: 320px;
    box-sizing: border-box;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    gap: 10px;
    background-color: #ffffff;
    padding: 14px 16px;
    border: 1px solid rgba(15, 23, 42, 0.07);
    border-radius: 12px;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}

.multi-modal-kde-container img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
}

.multi-modal-kde-label {
    font-size: 11px;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    color: #94a3b8;
    font-weight: 600;
    text-align: center;
    margin-bottom: 6px;
}

.multi-modal-img-container {
    max-width: 80%;
    width: 80%;
    min-width: 480px;
    height: 140px;
    padding: 14px 16px;
    border-radius: 12px;
    display: flex;
    justify-content: center;
    align-items: center;
    background-color: #ffffff;
    border: 1px solid rgba(15, 23, 42, 0.07);
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}

.multi-modal-img-container img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
    border-radius: 4px;
}

.multi-modal-slider {
    width: 480px;
    min-width: 480px;
    max-width: 480px;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 10px;
    padding: 14px 16px;
    border-radius: 12px;
    border: 1px solid rgba(15, 23, 42, 0.07);
    background-color: #ffffff;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
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
    box-sizing: border-box;
    padding: 10px 14px;
    margin: 6px 0;
    border-radius: 10px;
    font-size: 13px;
    line-height: 1.5;
    color: #334155;
    border: 1px solid rgba(15, 23, 42, 0.06);
    box-shadow: 0 1px 1px rgba(15, 23, 42, 0.03);
    transition: border-color 0.18s ease, box-shadow 0.18s ease;
}

.multi-info-item:hover {
    border-color: rgba(15, 23, 42, 0.12);
    box-shadow: 0 2px 6px rgba(15, 23, 42, 0.06);
}

.multi-info-item strong,
.multi-info-item b {
    color: #0f172a;
    font-weight: 600;
}

.multi-info-item a {
    color: #2563eb;
    text-decoration: none;
    font-weight: 500;
    padding: 2px 8px;
    border-radius: 999px;
    background-color: rgba(37, 99, 235, 0.08);
    transition: background-color 0.18s ease, color 0.18s ease;
}

.multi-info-item a:hover {
    background-color: rgba(37, 99, 235, 0.16);
    text-decoration: none;
}

.multi-info-item:empty {
    display: none;
}


"""
