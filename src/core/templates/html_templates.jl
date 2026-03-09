# ============================================================================
# HTML TEMPLATES
# ============================================================================

"""
HTML template for compound motifs (pairs, triplets, quadruplets).

Now supports multiple groups on the same page with custom button texts.

Mustache variables:
- `protein_name`: Title displayed in page
- `j`: Page number (for linking correct JS file)
- `DF`: DataFrame with columns:
  - div_img_id, i, img_src, img_alt
  - div_text_id, p_id1-6_default
  - div_slide_id, max_comb
  - group_id, button_text (NEW)

Features:
- Responsive grid layout
- Image/text containers for each motif
- Range sliders for distance selection
- Multiple collapsible groups with custom button labels
- Multiple modal windows (sequence highlight, cluster view, consensus, histograms)
"""
html_template = mt"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{{:protein_name}}} motifs</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>    
    <br><br>
    <div class="wrapper">

        <div id="nav" style="display: flex; justify-content: center;"></div>
        <br><br>
        
        <script>
        // Group multi-motif data by group_id and button_text
        const multiGroupedData = {};
        const multiFullData = [
            {{#:DF}}
            {
                div_img_id: "{{:div_img_id}}",
                i: {{:i}},
                img_src: "{{{:img_src}}}",
                img_alt: "{{:img_alt}}",
                div_text_id: "{{:div_text_id}}",
                p_id1: `{{{:p_id1_default}}}`,
                p_id2: `{{{:p_id2_default}}}`,
                p_id3: `{{{:p_id3_default}}}`,
                p_id4: `{{{:p_id4_default}}}`,
                p_id5: `{{{:p_id5_default}}}`,
                p_id6: `{{{:p_id6_default}}}`,
                p_id7: `{{{:p_id7_default}}}`,
                div_slide_id: "{{:div_slide_id}}",
                max_comb: {{:max_comb}},
                filter_indices: "{{{:filter_indices}}}",
                median: "{{:median_val}}",
                group_id: "{{:group_id}}",
                button_text: "{{:button_text}}",
                has_interaction: {{:has_interaction}}
            },
            {{/:DF}}
        ].filter(x => x.img_src);
        
        // Group by group_id
        multiFullData.forEach(item => {
            const gid = item.group_id || "default";
            if (!multiGroupedData[gid]) {
                multiGroupedData[gid] = {
                    button_text: item.button_text || "Multi-Motifs",
                    items: []
                };
            }
            multiGroupedData[gid].items.push(item);
        });
        </script>
        
        <div id="multiGroupsContainer"></div>
        
        <script>
        // Dynamically create toggle buttons and grids for each group
        const multiContainer = document.getElementById('multiGroupsContainer');
        let multiGroupIndex = 0;
        
        for (const [groupId, groupInfo] of Object.entries(multiGroupedData)) {
            multiGroupIndex++;
            const isFirstGroup = multiGroupIndex === 1;
            const wrapperId = 'multiGridWrapper_' + groupId;
            const buttonId = 'multiToggleButton_' + groupId;
            
            // Determine initial state based on group index
            const initiallyVisible = isFirstGroup;
            const buttonClass = initiallyVisible ? 'active' : '';
            const wrapperClass = initiallyVisible ? 'visible' : '';
            const iconText = initiallyVisible ? '▼' : '▶';
            const buttonTextContent = initiallyVisible ? groupInfo.button_text : 'Show ' + groupInfo.button_text;
            
            // Create toggle button section
            const toggleContainer = document.createElement('div');
            toggleContainer.className = 'grid-toggle-container';
            toggleContainer.innerHTML = `
                <div id="multiToggleBar_${groupId}" class="grid-toggle-bar ${buttonClass}">
                    <button id="${buttonId}" class="grid-toggle-button ${buttonClass}" onclick="toggleMultiGrid('${groupId}')">
                        <span class="grid-toggle-icon">${iconText}</span>
                        <span>${buttonTextContent}</span>
                    </button>
                </div>
            `;
            multiContainer.appendChild(toggleContainer);
            
            // Create grid wrapper
            const wrapper = document.createElement('div');
            wrapper.id = wrapperId;
            wrapper.className = 'grid-wrapper ' + wrapperClass;
            
            // Create container
            const container = document.createElement('div');
            container.className = 'container';
            
            // Add slider groups for each item in this group
            groupInfo.items.forEach(item => {
                const sliderGroup = document.createElement('div');
                sliderGroup.className = 'sliderGroup';
                sliderGroup.setAttribute('data-median', item.median);
                sliderGroup.onclick = () => openMultiModal(item.i);
                
                sliderGroup.innerHTML = `
                    <div class="imageTextContainer">
                        <div id="${item.div_img_id}" class="imageContainer">
                            <img id="img${item.i}" src="${item.img_src}" alt="${item.img_alt}">
                            <span class="filter-index-overlay">${item.filter_indices}</span>
                            ${item.has_interaction ? '<span class="interaction-indicator">I</span>' : ''}
                        </div>
                        <div id="${item.div_text_id}" class="textContainer">
                            ${item.p_id1 ? `<p id="text${item.i}_1" class="imageText">${item.p_id1}</p>` : ''}
                            ${item.p_id2 ? `<p id="text${item.i}_2" class="imageText">${item.p_id2}</p>` : ''}
                            ${item.p_id3 ? `<p id="text${item.i}_3" class="imageText">${item.p_id3}</p>` : ''}
                            ${item.p_id4 ? `<p id="text${item.i}_4" class="imageText">${item.p_id4}</p>` : ''}
                            ${item.p_id5 ? `<p id="text${item.i}_5" class="imageText">${item.p_id5}</p>` : ''}
                            ${item.p_id6 ? `<p id="text${item.i}_6" class="imageText">${item.p_id6}</p>` : ''}
                            ${item.p_id7 ? `<p id="text${item.i}_7" class="imageText">${item.p_id7}</p>` : ''}
                        </div>
                    </div>
                    <div id="${item.div_slide_id}" class="sliderContainer">
                        <input id="valR${item.i}" type="range" min="0" max="${item.max_comb}" value="0">
                        <span id="range${item.i}">Image 1</span>
                    </div>
                `;
                
                container.appendChild(sliderGroup);
            });
            
            wrapper.appendChild(container);
            multiContainer.appendChild(wrapper);
        }
        
        // Toggle function for multi-motif grids
        function toggleMultiGrid(groupId) {
            const wrapper = document.getElementById('multiGridWrapper_' + groupId);
            const button = document.getElementById('multiToggleButton_' + groupId);
            const buttonText = button.querySelector('span:last-child');
            const icon = button.querySelector('.grid-toggle-icon');
            
            // Find original button text from groupedData
            const originalText = multiGroupedData[groupId].button_text;
            
            if (wrapper.classList.contains('visible')) {
                wrapper.classList.remove('visible');
                button.classList.remove('active');
                buttonText.textContent = 'Show ' + originalText;
                icon.textContent = '▶';
            } else {
                wrapper.classList.add('visible');
                button.classList.add('active');
                buttonText.textContent = originalText;
                icon.textContent = '▼';
            }
        }
        
        // Populate multiModalData from jsonData for modal functionality
        // This needs to happen after data{{:j}}.js is loaded
        window.addEventListener('DOMContentLoaded', function() {
            // Apply dynamic card styling for colored borders
            if (typeof applyDynamicCardStyling === 'function') {
                applyDynamicCardStyling();
            }
            
            if (typeof jsonData !== 'undefined') {
                // Build arrays in index order (sorted by i)
                const sortedData = multiFullData.sort((a, b) => a.i - b.i);
                const images = [];
                const labels = [];
                const texts = [];
                const baseFolders = [];
                
                sortedData.forEach(item => {
                    // Get mode key for this item
                    const modeKey = 'mode_' + (item.group_id ? item.group_id + '_' : '') + item.i;
                    
                    if (jsonData[modeKey]) {
                        images.push(jsonData[modeKey].pwms || []);
                        labels.push(jsonData[modeKey].labels || []);
                        texts.push(jsonData[modeKey].texts || []);
                        
                        // Extract base folder from first image
                        const pwms = jsonData[modeKey].pwms;
                        if (pwms && pwms.length > 0) {
                            const firstImage = pwms[0];
                            const lastSlashIndex = firstImage.lastIndexOf('/');
                            baseFolders.push(lastSlashIndex >= 0 ? firstImage.substring(0, lastSlashIndex) : '');
                        } else {
                            baseFolders.push('');
                        }
                    }
                });
                
                // Store for modal access
                if (typeof multiModalData !== 'undefined') {
                    multiModalData.images = images;
                    multiModalData.labels = labels;
                    multiModalData.texts = texts;
                    multiModalData.baseFolders = baseFolders;
                }
            }
        });
        </script>
    </div>
    <script src="data{{:j}}.js"></script>
    <script src="scripts{{:j}}.js"></script>
    <!-- Modal Structure -->
    <div id="highlightModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal()">&times;</span>
            <div id="highlightedSequences"></div>
        </div>
    </div>

    <!-- Modal Structure for clustering -->
    <div id="highlightModal_cluster">
        <div id="highlightContent">
            <span class="close" onclick="closeModal_cluster()">&times;</span>
            <div class="modal-column">
                <img id="modalImage" src="" alt="Image">
            </div>
            <div class="modal-column">
                <p id="modalText">This is placeholder text. You can replace it with dynamic content as needed.</p>
            </div>
        </div>
    </div>

    <!-- Modal Structure for consensus_str -->
    <div id="highlightModal_text">
        <div id="highlightModal_text_content">
            <span class="close" onclick="closeModal_text()">&times;</span>
            <div class="modal-column">
                <p id="modalText1"> This is placeholder text. You can replace it with dynamic content as needed.</p>
                <button id="copyButton" onclick="copyText()">copy string</button>
            </div>
        </div>
    </div>

    <!-- Modal Structure for just image -->
    <div id="highlightModal_img">
        <div id="highlightModal_img_content">
            <span class="close" onclick="closeModal_img()">&times;</span>
            <div class="modal-column">
                <img id="modalImage1" src="" alt="Image">
            </div>
        </div>
    </div>

    <!-- Modal Structure for multi-motif view -->
    <div id="multiMotifModal" class="multi-modal">
        <div class="multi-modal-content">
            <span class="multi-close" onclick="closeMultiMotifModal()">&times;</span>
            <div class="multi-modal-body">
                <div class="multi-modal-left">
                    <div class="multi-modal-influence-container">
                        <span class="multi-modal-influence-label">Fixed Distance Influence</span>
                        <img id="multiMotifInfluenceFixed" src="" alt="Fixed Distance Influence">
                    </div>
                    <div class="multi-modal-influence-container">
                        <span class="multi-modal-influence-label">Relaxed Distance Influence</span>
                        <img id="multiMotifInfluenceRelaxed" src="" alt="Relaxed Distance Influence">
                    </div>
                    <div class="multi-modal-info">
                        <div id="multiMotifText1" class="multi-info-item"></div>
                        <div id="multiMotifText2" class="multi-info-item"></div>
                        <div id="multiMotifText3" class="multi-info-item"></div>
                        <div id="multiMotifText4" class="multi-info-item"></div>
                        <div id="multiMotifText5" class="multi-info-item"></div>
                        <div id="multiMotifText6" class="multi-info-item"></div>
                    </div>
                </div>
                <div class="multi-modal-right">
                    <div class="multi-modal-kde-container">
                        <span class="multi-modal-kde-label">Indicator Plot</span>
                        <img id="multiMotifKde" src="" alt="KDE Plot">
                    </div>
                    <div class="multi-modal-img-container">
                        <img id="multiMotifImage" src="" alt="Motif Image">
                    </div>
                    <div class="multi-modal-slider">
                        <input id="multiMotifSlider" type="range" min="0" max="100" value="0">
                        <span id="multiMotifRangeLabel" class="multi-modal-slider-label">Position: 0</span>
                    </div>
                </div>
            </div>
        </div>
    </div>
"""

"""
HTML template for singleton motifs.

Similar to html_template but without sliders (singletons have no distance variations).
Now supports multiple groups on the same page with custom button texts.

Mustache variables:
- `protein_name`: Title displayed in page
- `j`: Page number
- `DF`: DataFrame with motif data (includes group_id and button_text columns)

Features:
- Simpler layout (no sliders)
- Static image display
- Hover-enlarge and click-to-expand modal
- Multiple collapsible groups with custom button labels
"""
html_template_singleton = mt"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{{:protein_name}}} motifs</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <br><br>
    <div class="wrapper">
        <div id="nav" style="display: flex; justify-content: center;"></div>
        <br><br>
        
        <script>
        // Group singleton data by group_id and button_text
        const groupedData = {};
        const singletonFullData = [
            {{#:DF}}
            {
                img: "{{{:img_src}}}",
                influence: "{{{:img_src}}}".replace('.png', '_influence.png'),
                yy_kde: "{{{:img_src}}}".replace(/\/([^\/]+)\.png$/, '/yy_kde_intersect_$1.png'),
                alt: "{{:img_alt}}",
                filter_indices: "{{{:filter_indices}}}",
                text1: `{{{:p_id1_default}}}`,
                text2: `{{{:p_id2_default}}}`,
                text3: `{{{:p_id3_default}}}`,
                text4: `{{{:p_id4_default}}}`,
                text5: `{{{:p_id5_default}}}`,
                text6: `{{{:p_id6_default}}}`,
                text7: `{{{:p_id7_default}}}`,
                median: "{{:median_val}}",
                group_id: "{{:group_id}}",
                button_text: "{{:button_text}}",
                i: {{:i}},
                is_significant: {{:is_significant}}
            },
            {{/:DF}}
        ].filter(x => x.img);
        
        // Group by group_id
        singletonFullData.forEach(item => {
            const gid = item.group_id || "default";
            if (!groupedData[gid]) {
                groupedData[gid] = {
                    button_text: item.button_text || "Singleton Motifs",
                    items: []
                };
            }
            groupedData[gid].items.push(item);
        });
        
        // Store for modal access
        window.singletonData = singletonFullData;
        </script>
        
        <div id="singletonGroupsContainer"></div>
        
        <script>
        // Dynamically create toggle buttons and grids for each group
        const container = document.getElementById('singletonGroupsContainer');
        let groupIndex = 0;
        
        for (const [groupId, groupInfo] of Object.entries(groupedData)) {
            groupIndex++;
            const isFirstGroup = groupIndex === 1;
            const wrapperId = 'singletonGridWrapper_' + groupId;
            const buttonId = 'singletonToggleButton_' + groupId;
            
            // Determine initial state based on group index
            const initiallyVisible = isFirstGroup;
            const buttonClass = initiallyVisible ? 'active' : '';
            const wrapperClass = initiallyVisible ? 'visible' : '';
            const iconText = initiallyVisible ? '▼' : '▶';
            const buttonTextContent = initiallyVisible ? groupInfo.button_text : 'Show ' + groupInfo.button_text;
            
            // Create toggle button section
            const toggleContainer = document.createElement('div');
            toggleContainer.className = 'grid-toggle-container';
            toggleContainer.innerHTML = `
                <div id="singletonToggleBar_${groupId}" class="grid-toggle-bar ${buttonClass}">
                    <button id="${buttonId}" class="grid-toggle-button ${buttonClass}" onclick="toggleSingletonGrid('${groupId}')">
                        <span class="grid-toggle-icon">${iconText}</span>
                        <span>${buttonTextContent}</span>
                    </button>
                </div>
            `;
            container.appendChild(toggleContainer);
            
            // Create grid wrapper
            const wrapper = document.createElement('div');
            wrapper.id = wrapperId;
            wrapper.className = 'grid-wrapper ' + wrapperClass;
            
            // Create grid
            const grid = document.createElement('div');
            grid.className = 'singleton-grid';
            
            // Add cells for each item in this group
            groupInfo.items.forEach(item => {
                const cell = document.createElement('div');
                cell.className = 'singleton-cell';
                cell.setAttribute('data-median', item.median);
                cell.onclick = () => openSingletonModal(item.i);
                
                cell.innerHTML = `
                    <img src="${item.img}" alt="${item.alt}" class="singleton-img">
                    <span class="singleton-filter-overlay">${item.alt.replace('motif ', '')}</span>
                    ${!item.is_significant ? '<span class="insignificant-indicator">insignificant</span>' : ''}
                `;
                
                grid.appendChild(cell);
            });
            
            wrapper.appendChild(grid);
            container.appendChild(wrapper);
        }
        
        // Toggle function for singleton grids
        function toggleSingletonGrid(groupId) {
            const wrapper = document.getElementById('singletonGridWrapper_' + groupId);
            const button = document.getElementById('singletonToggleButton_' + groupId);
            const buttonText = button.querySelector('span:last-child');
            const icon = button.querySelector('.grid-toggle-icon');
            
            // Find original button text from groupedData
            const originalText = groupedData[groupId].button_text;
            
            if (wrapper.classList.contains('visible')) {
                wrapper.classList.remove('visible');
                button.classList.remove('active');
                buttonText.textContent = 'Show ' + originalText;
                icon.textContent = '▶';
            } else {
                wrapper.classList.add('visible');
                button.classList.add('active');
                buttonText.textContent = originalText;
                icon.textContent = '▼';
            }
        }
        
        // Apply dynamic card styling when DOM is ready
        window.addEventListener('DOMContentLoaded', function() {
            if (typeof applyDynamicCardStyling === 'function') {
                applyDynamicCardStyling();
            }
        });
        </script>
    </div>
    
    <!-- Singleton Modal -->
    <div id="singletonModal" class="singleton-modal">
        <div class="singleton-modal-content">
            <span class="singleton-close" onclick="closeSingletonModal()">&times;</span>
            <div class="singleton-modal-body">
                <div class="singleton-modal-left">
                    <div class="singleton-modal-influence-container">
                        <span class="singleton-modal-influence-label">Influence</span>
                        <img id="singletonModalInfluence" src="" alt="Influence Plot">
                    </div>
                    <div class="singleton-modal-info">
                        <div id="singletonModalText1" class="singleton-info-item"></div>
                        <div id="singletonModalText2" class="singleton-info-item"></div>
                        <div id="singletonModalText3" class="singleton-info-item"></div>
                        <div id="singletonModalText4" class="singleton-info-item"></div>
                        <div id="singletonModalText5" class="singleton-info-item"></div>
                        <div id="singletonModalText6" class="singleton-info-item"></div>
                    </div>
                </div>
                <div class="singleton-modal-right">
                    <div class="singleton-modal-kde-container">
                        <span class="singleton-modal-kde-label">Indicator Plot</span>
                        <img id="singletonModalKde" src="" alt="KDE Plot">
                    </div>
                    <div class="singleton-modal-img-container">
                        <img id="singletonModalImg" src="" alt="">
                    </div>
                    <h3 id="singletonModalTitle"></h3>
                </div>
            </div>
        </div>
    </div>

    <script src="data{{:j}}.js"></script>
    <script src="scripts{{:j}}.js"></script>

    <!-- Modal Structure for sequence highlighting -->
    <div id="highlightModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal()">&times;</span>
            <div id="highlightedSequences"></div>
        </div>
    </div>
"""

"""
HTML template combining both singletons and multi-motifs on the same page.

This unified template allows you to mix singleton and multi-motif groups with
custom button labels for each group.

Mustache variables:
- `protein_name`: Title displayed in page
- `j`: Page number
- `DF`: DataFrame with ALL motif data (both singletons and multi-motifs)
  Must include: group_id, button_text, max_comb (0 for singletons, >0 for multi)

Features:
- Combines singleton and multi-motif displays
- Multiple collapsible groups with custom labels
- Automatic detection of singleton vs multi-motif based on max_comb value
- All modals and interactions work seamlessly
"""
html_template_unified = mt"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{{:protein_name}}} motifs</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <br><br>
    <div class="wrapper">
        <div id="nav" style="display: flex; justify-content: center;"></div>
        <br><br>
        
        <script>
        // All motif data (both singletons and multi-motifs)
        const allMotifData = [
            {{#:DF}}
            {
                // Common fields
                i: {{:i}},
                img_src: "{{{:img_src}}}",
                img: "{{{:img_src}}}",  // Alias for singleton modal compatibility
                img_alt: "{{:img_alt}}",
                alt: "{{:img_alt}}",  // Alias for singleton modal compatibility
                text1: `{{{:p_id1_default}}}`,
                text2: `{{{:p_id2_default}}}`,
                text3: `{{{:p_id3_default}}}`,
                text4: `{{{:p_id4_default}}}`,
                text5: `{{{:p_id5_default}}}`,
                text6: `{{{:p_id6_default}}}`,
                text7: `{{{:p_id7_default}}}`,
                median: "{{:median_val}}",
                group_id: "{{:group_id}}",
                button_text: "{{:button_text}}",
                filter_indices: "{{{:filter_indices}}}",
                // Multi-motif specific
                max_comb: {{:max_comb}},
                div_img_id: "{{:div_img_id}}",
                div_text_id: "{{:div_text_id}}",
                div_slide_id: "{{:div_slide_id}}",
                // Derived
                influence: "{{{:img_src}}}".replace('.png', '_influence.png'),
                yy_kde: "{{{:img_src}}}".replace(/\/([^\/]+)\.png$/, '/yy_kde_intersect_$1.png'),
                is_singleton: {{:max_comb}} === 0,
                has_interaction: {{:has_interaction}},
                is_significant: {{:is_significant}}
            },
            {{/:DF}}
        ].filter(x => x.img_src);
        
        // Group by group_id
        const motifGroupedData = {};
        allMotifData.forEach(item => {
            const gid = item.group_id || "default";
            if (!motifGroupedData[gid]) {
                motifGroupedData[gid] = {
                    button_text: item.button_text || "Motifs",
                    is_singleton: item.is_singleton,
                    items: []
                };
            }
            motifGroupedData[gid].items.push(item);
        });
        
        // Store for modal access
        window.singletonData = allMotifData.filter(x => x.is_singleton);
        </script>
        
        <div id="unifiedGroupsContainer"></div>
        
        <script>
        // Dynamically create groups
        const unifiedContainer = document.getElementById('unifiedGroupsContainer');
        
        // Track group index for default visibility (first group visible, rest hidden)
        let groupIndex = 0;
        
        for (const [groupId, groupInfo] of Object.entries(motifGroupedData)) {
            groupIndex++;
            const isFirstGroup = groupIndex === 1;
            const wrapperId = 'gridWrapper_' + groupId;
            const buttonId = 'toggleButton_' + groupId;
            const isSingleton = groupInfo.is_singleton;
            
            // Determine initial state based on group index
            const initiallyVisible = isFirstGroup;
            const buttonClass = initiallyVisible ? 'active' : '';
            const wrapperClass = initiallyVisible ? 'visible' : '';
            const iconText = initiallyVisible ? '▼' : '▶';
            const buttonTextContent = initiallyVisible ? groupInfo.button_text : 'Show ' + groupInfo.button_text;
            
            // Create toggle button
            const toggleContainer = document.createElement('div');
            toggleContainer.className = 'grid-toggle-container';
            toggleContainer.innerHTML = `
                <div class="grid-toggle-bar ${buttonClass}">
                    <button id="${buttonId}" class="grid-toggle-button ${buttonClass}" onclick="toggleGrid('${groupId}')">
                        <span class="grid-toggle-icon">${iconText}</span>
                        <span>${buttonTextContent}</span>
                    </button>
                </div>
            `;
            unifiedContainer.appendChild(toggleContainer);
            
            // Create grid wrapper
            const wrapper = document.createElement('div');
            wrapper.id = wrapperId;
            wrapper.className = 'grid-wrapper ' + wrapperClass;
            
            if (isSingleton) {
                // Singleton grid
                const grid = document.createElement('div');
                grid.className = 'singleton-grid';
                groupInfo.items.forEach(item => {
                    const cell = document.createElement('div');
                    cell.className = 'singleton-cell';
                    cell.setAttribute('data-median', item.median);
                    cell.onclick = () => openSingletonModal(item.i);
                    cell.innerHTML = `
                        <img src="${item.img_src}" alt="${item.img_alt}" class="singleton-img">
                        <span class="singleton-filter-overlay">${item.filter_indices}</span>
                        ${!item.is_significant ? '<span class="insignificant-indicator">insignificant</span>' : ''}
                    `;
                    grid.appendChild(cell);
                });
                wrapper.appendChild(grid);
            } else {
                // Multi-motif grid
                const container = document.createElement('div');
                container.className = 'container';
                groupInfo.items.forEach(item => {
                    const sliderGroup = document.createElement('div');
                    sliderGroup.className = 'sliderGroup';
                    sliderGroup.setAttribute('data-median', item.median);
                    sliderGroup.onclick = () => openMultiModal(item.i);
                    sliderGroup.innerHTML = `
                        <div class="imageTextContainer">
                            <div id="${item.div_img_id}" class="imageContainer">
                                <img id="img${item.i}" src="${item.img_src}" alt="${item.img_alt}">
                                <span class="filter-index-overlay">${item.filter_indices}</span>
                                ${item.has_interaction ? '<span class="interaction-indicator">I</span>' : ''}
                            </div>
                            <div id="${item.div_text_id}" class="textContainer">
                                ${item.text1 ? `<p id="text${item.i}_1" class="imageText">${item.text1}</p>` : ''}
                                ${item.text2 ? `<p id="text${item.i}_2" class="imageText">${item.text2}</p>` : ''}
                                ${item.text3 ? `<p id="text${item.i}_3" class="imageText">${item.text3}</p>` : ''}
                                ${item.text4 ? `<p id="text${item.i}_4" class="imageText">${item.text4}</p>` : ''}
                                ${item.text5 ? `<p id="text${item.i}_5" class="imageText">${item.text5}</p>` : ''}
                                ${item.text6 ? `<p id="text${item.i}_6" class="imageText">${item.text6}</p>` : ''}
                                ${item.text7 ? `<p id="text${item.i}_7" class="imageText">${item.text7}</p>` : ''}
                            </div>
                        </div>
                        <div id="${item.div_slide_id}" class="sliderContainer">
                            <input id="valR${item.i}" type="range" min="0" max="${item.max_comb}" value="0">
                            <span id="range${item.i}">Image 1</span>
                        </div>
                    `;
                    container.appendChild(sliderGroup);
                });
                wrapper.appendChild(container);
            }
            
            unifiedContainer.appendChild(wrapper);
        }
        
        // Toggle function
        function toggleGrid(groupId) {
            const wrapper = document.getElementById('gridWrapper_' + groupId);
            const button = document.getElementById('toggleButton_' + groupId);
            const buttonText = button.querySelector('span:last-child');
            const icon = button.querySelector('.grid-toggle-icon');
            const originalText = motifGroupedData[groupId].button_text;
            
            if (wrapper.classList.contains('visible')) {
                wrapper.classList.remove('visible');
                button.classList.remove('active');
                buttonText.textContent = 'Show ' + originalText;
                icon.textContent = '▶';
            } else {
                wrapper.classList.add('visible');
                button.classList.add('active');
                buttonText.textContent = originalText;
                icon.textContent = '▼';
            }
        }
        
        // Populate multiModalData from jsonData for modal functionality
        // This needs to happen after data{{:j}}.js is loaded
        window.addEventListener('DOMContentLoaded', function() {
            // Apply dynamic card styling for colored borders
            if (typeof applyDynamicCardStyling === 'function') {
                applyDynamicCardStyling();
            }
            
            if (typeof jsonData !== 'undefined') {
                // Build arrays in index order (sorted by i) for multi-motifs only
                const multiMotifData = allMotifData.filter(x => !x.is_singleton).sort((a, b) => a.i - b.i);
                const images = [];
                const labels = [];
                const texts = [];
                const baseFolders = [];
                
                // Create mapping from display index (i) to array index
                window.multiMotifIndexMap = {};
                
                multiMotifData.forEach((item, arrayIndex) => {
                    // Map display index to array index
                    window.multiMotifIndexMap[item.i] = arrayIndex;
                    
                    // Get mode key for this item
                    const modeKey = 'mode_' + (item.group_id ? item.group_id + '_' : '') + item.i;
                    
                    if (jsonData[modeKey]) {
                        images.push(jsonData[modeKey].pwms || []);
                        labels.push(jsonData[modeKey].labels || []);
                        texts.push(jsonData[modeKey].texts || []);
                        
                        // Extract base folder from first image
                        const pwms = jsonData[modeKey].pwms;
                        if (pwms && pwms.length > 0) {
                            const firstImage = pwms[0];
                            const lastSlashIndex = firstImage.lastIndexOf('/');
                            baseFolders.push(lastSlashIndex >= 0 ? firstImage.substring(0, lastSlashIndex) : '');
                        } else {
                            baseFolders.push('');
                        }
                    }
                });
                
                // Store for modal access
                if (typeof multiModalData !== 'undefined') {
                    multiModalData.images = images;
                    multiModalData.labels = labels;
                    multiModalData.texts = texts;
                    multiModalData.baseFolders = baseFolders;
                }
            }
        });
        </script>
    </div>
    
    <!-- All modals -->
    <div id="singletonModal" class="singleton-modal">
        <div class="singleton-modal-content">
            <span class="singleton-close" onclick="closeSingletonModal()">&times;</span>
            <div class="singleton-modal-body">
                <div class="singleton-modal-left">
                    <div class="singleton-modal-influence-container">
                        <span class="singleton-modal-influence-label">Influence</span>
                        <img id="singletonModalInfluence" src="" alt="Influence Plot">
                    </div>
                    <div class="singleton-modal-info">
                        <div id="singletonModalText1" class="singleton-info-item"></div>
                        <div id="singletonModalText2" class="singleton-info-item"></div>
                        <div id="singletonModalText3" class="singleton-info-item"></div>
                        <div id="singletonModalText4" class="singleton-info-item"></div>
                        <div id="singletonModalText5" class="singleton-info-item"></div>
                        <div id="singletonModalText6" class="singleton-info-item"></div>
                    </div>
                </div>
                <div class="singleton-modal-right">
                    <div class="singleton-modal-kde-container">
                        <span class="singleton-modal-kde-label">Indicator Plot</span>
                        <img id="singletonModalKde" src="" alt="KDE Plot">
                    </div>
                    <div class="singleton-modal-img-container">
                        <img id="singletonModalImg" src="" alt="">
                    </div>
                    <h3 id="singletonModalTitle"></h3>
                </div>
            </div>
        </div>
    </div>
    
    <div id="multiMotifModal" class="multi-modal">
        <div class="multi-modal-content">
            <span class="multi-close" onclick="closeMultiMotifModal()">&times;</span>
            <div class="multi-modal-body">
                <div class="multi-modal-left">
                    <div class="multi-modal-influence-container">
                        <span class="multi-modal-influence-label">Fixed Distance Influence</span>
                        <img id="multiMotifInfluenceFixed" src="" alt="Fixed Distance Influence">
                    </div>
                    <div class="multi-modal-influence-container">
                        <span class="multi-modal-influence-label">Relaxed Distance Influence</span>
                        <img id="multiMotifInfluenceRelaxed" src="" alt="Relaxed Distance Influence">
                    </div>
                    <div class="multi-modal-info">
                        <div id="multiMotifText1" class="multi-info-item"></div>
                        <div id="multiMotifText2" class="multi-info-item"></div>
                        <div id="multiMotifText3" class="multi-info-item"></div>
                        <div id="multiMotifText4" class="multi-info-item"></div>
                        <div id="multiMotifText5" class="multi-info-item"></div>
                        <div id="multiMotifText6" class="multi-info-item"></div>
                    </div>
                </div>
                <div class="multi-modal-right">
                    <div class="multi-modal-kde-container">
                        <span class="multi-modal-kde-label">Indicator Plot</span>
                        <img id="multiMotifKde" src="" alt="KDE Plot">
                    </div>
                    <div class="multi-modal-img-container">
                        <img id="multiMotifImage" src="" alt="Motif Image">
                    </div>
                    <div class="multi-modal-slider">
                        <input id="multiMotifSlider" type="range" min="0" max="100" value="0">
                        <span id="multiMotifRangeLabel" class="multi-modal-slider-label">Position: 0</span>
                    </div>
                </div>
            </div>
        </div>
    </div>
    
    <div id="highlightModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal()">&times;</span>
            <div id="highlightedSequences"></div>
        </div>
    </div>

    <script src="data{{:j}}.js"></script>
    <script src="scripts{{:j}}.js"></script>
"""

"""
Hover window template for displaying metadata.

Mustache variables:
- `hover_on`: HTML attributes for hover element
- `meta_str`: Metadata content to display
"""
html_hover_default = mt"""
    <div {{{:hover_on}}}>
        <div class="hover-meta-data">
        {{{:meta_str}}}
        </div>
    </div>"""

"""
HTML template for the generalization page (index2.html).

Displays the generalization scatter plot image centered on the page
with the same navigation bar as the main motif pages.

Mustache variables:
- `protein_name`: Title displayed in page
- `j`: Page number (should be 2 for generalization)
- `upto`: Total number of navigation pages
- `image_path`: Path to the generalization.png image
"""
html_template_generalization = mt"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{{:protein_name}}} - Generalization</title>
    <link rel="stylesheet" href="styles.css">
    <style>
        .generalization-container {
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            min-height: 60vh;
            padding: 40px 20px;
        }
        .generalization-title {
            font-size: 24px;
            font-weight: 600;
            color: #333;
            margin-bottom: 30px;
            text-align: center;
            font-family: 'Helvetica', 'Arial', sans-serif;
        }
        .generalization-image {
            max-width: 95%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        .generalization-caption {
            margin-top: 20px;
            font-size: 14px;
            color: #666;
            text-align: center;
            max-width: 800px;
        }
    </style>
</head>
<body>
    <br><br>
    <div class="wrapper">
        <div id="nav" style="display: flex; justify-content: center;"></div>
        <br><br>
        
        <div class="generalization-container">
            <!--- <h1 class="generalization-title">Generalization on the test set.</h1> --->
            <img src="{{{:image_path}}}" alt="Generalization Scatter Plots" class="generalization-image">
            <p class="generalization-caption">
                Generalization on the test set. Scatter plots comparing model predictions, labels, and processor products. 
                R² values indicate the correlation strength between each pair.
            </p>
        </div>
    </div>
    
    <script>
    function createArray(num) {
        return Array.from({ length: num }, (_, i) => i + 1);
    }
    
    function updateNav() {
        const currentPage = 'index{{:j}}.html';
        const availablePages = createArray({{:upto}});
        
        const navHtml = availablePages.map(num => {
            const page = 'index' + num + '.html';
            const label = num === 1 ? 'Motif influence' :
                          num === 2 ? 'Generalization' :
                          num === 3 ? 'Statistics' :
                          num === 4 ? 'Readme' :
                          num === 5 ? 'Page 5' :
                          'Page 6';
            
            return '<a href="' + page + '" ' + (currentPage === page ? 'class="current"' : '') + '>' + label + '</a>';
        }).join(' &nbsp;&nbsp; | &nbsp;&nbsp; ');
        
        document.getElementById('nav').innerHTML = navHtml;
    }
    
    window.onload = updateNav;
    </script>
</body>
</html>
"""

"""Closing HTML tags"""
html_end = "</body></html>"
