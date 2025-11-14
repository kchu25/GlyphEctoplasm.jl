# ============================================================================
# JAVASCRIPT STYLING FUNCTIONS
# ============================================================================

"""
JavaScript functions for dynamic card styling based on Banzhaf values.

Functions:
- `applyDynamicCardStyling()`: Main function to apply styling to all cards
- `styleCardBorders()`: Apply gradient border effects based on median values

Features:
1. Border thickness varies with absolute value of median (larger influence = thicker border)
2. Border color based on median value:
   - Positive values: dark red gradient
   - Negative values: dark blue gradient
   - Near-zero values: minimal border
3. Gradient effect: multiple layered box-shadows create outward fade
   - Inner border is most opaque
   - Outer layers become progressively more transparent
   - Creates a glowing halo effect that doesn't obscure adjacent cards
"""
script_styling_str = raw"""
// ============================================================================
// DYNAMIC CARD STYLING BASED ON MEDIAN VALUES
// ============================================================================

/**
 * Apply dynamic border styling to cards based on their median Banzhaf values.
 */
function applyDynamicCardStyling() {
    // Check if colored borders are enabled (can be disabled via global flag)
    if (typeof window.ENABLE_COLORED_BORDERS !== 'undefined' && !window.ENABLE_COLORED_BORDERS) {
        return;  // Skip styling if disabled
    }
    
    // Style multi-motif cards
    const multiCards = document.querySelectorAll('.sliderGroup[data-median]');
    const multiMedians = Array.from(multiCards).map(card => parseFloat(card.dataset.median));
    
    if (multiMedians.length > 0) {
        styleCardBorders(multiCards, multiMedians);
    }
    
    // Style singleton cards
    const singletonCards = document.querySelectorAll('.singleton-cell[data-median]');
    const singletonMedians = Array.from(singletonCards).map(card => parseFloat(card.dataset.median));
    
    if (singletonMedians.length > 0) {
        styleCardBorders(singletonCards, singletonMedians);
    }
}

/**
 * Apply gradient border styling to a set of cards.
 * 
 * @param {NodeList} cards - DOM elements to style
 * @param {Array} medians - Array of median values corresponding to cards
 */
function styleCardBorders(cards, medians) {
    // Find absolute max for normalization
    const absMedians = medians.map(Math.abs);
    const maxAbsMedian = Math.max(...absMedians);
    
    if (maxAbsMedian === 0) return;  // No styling needed if all are zero
    
    // Border parameters
    const minBorderWidth = 1;      // Minimum border width in pixels
    const maxBorderWidth = 4;     // Maximum border width in pixels (increased for emphasis)
    const numGradientLayers = 48;   // Number of gradient layers for smooth fade
    
    cards.forEach((card, index) => {
        const median = medians[index];
        const absMedian = Math.abs(median);
        
        // Use power function to emphasize larger values
        const normalizedAbs = Math.pow(absMedian / maxAbsMedian, 0.5);  // Power < 1 for more contrast
        
        // Calculate border width
        const borderWidth = minBorderWidth + normalizedAbs * (maxBorderWidth - minBorderWidth);
        
        // Determine base color based on sign of median
        let baseColor;
        if (Math.abs(median) < 0.01) {
            // Near-zero: very light gray
            baseColor = { r: 200, g: 200, b: 200 };
        } else if (median > 0) {
            // Positive: dark red
            baseColor = { r: 139, g: 0, b: 0 };
        } else {
            // Negative: Dark Blue
            baseColor = { r: 0, g: 0, b: 139 };
        }
        
        // Build layered box-shadow for gradient effect
        const boxShadows = [];
        
        for (let i = 0; i < numGradientLayers; i++) {
            // Calculate spread for this layer (outward from card)
            const layerSpread = (borderWidth / numGradientLayers) * (i + 1);
            
            // Calculate opacity that fades from inner to outer
            const maxOpacity = 0.05;
            const minOpacity = 0.01;
            const opacityFactor = 1 - (i / numGradientLayers);  // 1.0 -> 0.0
            const opacity = minOpacity + opacityFactor * (maxOpacity - minOpacity);
            
            // Scale opacity by normalized value
            const scaledOpacity = opacity * normalizedAbs;
            
            // Create rgba color string
            const color = `rgba(${baseColor.r}, ${baseColor.g}, ${baseColor.b}, ${scaledOpacity})`;
            
            // Add this layer to the shadow stack
            boxShadows.push(`0 0 ${Math.round(layerSpread * 2)}px ${Math.round(layerSpread)}px ${color}`);
        }
        
        // Apply solid border with base color
        const borderOpacity = 0.95 * normalizedAbs;
        const borderColor = `rgba(${baseColor.r}, ${baseColor.g}, ${baseColor.b}, ${borderOpacity})`;
        card.style.border = `${Math.max(1, Math.round(borderWidth / 3))}px solid ${borderColor}`;
        
        // Apply all gradient layers as box-shadow
        card.style.boxShadow = boxShadows.join(', ');
    });
}
"""
