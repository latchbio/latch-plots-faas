from agents import function_tool

@function_tool
async def export_to_csv(cell_id: str, filename: str) -> str:
    """Export cell output to CSV file.

    Args:
        cell_id: The cell containing dataframe output
        filename: Output filename (e.g., "results.csv")

    Returns:
        Success message with file path
    """
    return f"Exported to {filename}"

custom_tools = [
    export_to_csv,
]
