
# ============================================================================
# NUCLEOTIDE ALPHABET MAPPINGS
# ============================================================================

"""Mapping from nucleotide index to DNA character"""
const _ind2dna_str_ = Dict{Int, Char}(
    1 => 'A',
    2 => 'C',
    3 => 'G',
    4 => 'T'
)

"""Mapping from nucleotide index to RNA character"""
const _ind2dna_str_rna = Dict{Int, Char}(
    1 => 'A',
    2 => 'C',
    3 => 'G',
    4 => 'U'
)


"""
Placeholder character for low-confidence positions in consensus sequences.
Displayed when probability is below threshold.
"""
const _placeholder_char_ = 'n'


# ============================================================================
# JSON/HTML DATA STRUCTURE KEYS
# ============================================================================

"""JSON field name for PWM image paths"""
const pwms_str = "pwms"

"""JSON field name for descriptive labels"""
const labels_str = "labels"

"""JSON field name for text content arrays"""
const texts_str = "texts"

# ============================================================================
# HTML TAG IDENTIFIERS
# ============================================================================

"""HTML tag for image container div ID"""
const tag_div_img_id = "div_img_id"

"""HTML tag for index/iteration number"""
const tag_i = "i"

"""HTML tag for image source path"""
const tag_img_src = "img_src"

"""HTML tag for image alt text"""
const tag_img_alt = "img_alt"

"""HTML tag for text container div ID"""
const tag_div_text_id = "div_text_id"

"""HTML tag for first paragraph (contribution)"""
const tag_p_id1_default = "p_id1_default"

"""HTML tag for second paragraph (occurrences)"""
const tag_p_id2_default = "p_id2_default"

"""HTML tag for third paragraph (MEME link)"""
const tag_p_id3_default = "p_id3_default"

"""HTML tag for fourth paragraph (CSV link)"""
const tag_p_id4_default = "p_id4_default"

"""HTML tag for fifth paragraph (consensus)"""
const tag_p_id5_default = "p_id5_default"

"""HTML tag for sixth paragraph (gap histogram)"""
const tag_p_id6_default = "p_id6_default"

"""HTML tag for slider container div ID"""
const tag_div_slide_id = "div_slide_id"

"""HTML tag for maximum combination count"""
const tag_max_comb = "max_comb"

"""HTML tag for filter indices overlay"""
const tag_filter_indices = "filter_indices"

"""HTML tag for median Banzhaf value"""
const tag_median_val = "median_val"

"""HTML tag for group ID (for grouping singleton motifs)"""
const tag_group_id = "group_id"

"""HTML tag for button text (custom toggle button label)"""
const tag_button_text = "button_text"
