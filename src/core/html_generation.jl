"""
    fill_csv_link(csv_link) -> String

Create HTML link that opens sequence highlighting modal.

# Arguments
- `csv_link`: Path to CSV file containing sequences

# Returns
HTML anchor tag with onclick handler
"""
fill_csv_link(csv_link) =
    "<a href=\"#\" onclick=\"openHighlightPage(\'" * csv_link * "\')\">string highlight</a>"