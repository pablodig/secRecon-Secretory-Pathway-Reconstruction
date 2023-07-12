def get_gene_color(gene, gene_dict, process_dict, category_colors):

    """
    Function to determine the color of a gene based on its associated processes.

    It takes into consideration all processes associated with the gene, and assigns 
    a color based on the parent process that has the majority of processes associated 
    with the gene.

    Args:
    gene (str): The gene for which the color needs to be determined.
    gene_dict (dict): A dictionary with genes as keys and associated processes as values.
    process_dict (dict): A dictionary with categories (parent processes) as keys and associated processes as values.
    category_colors (dict): A dictionary with categories as keys and associated colors as values.

    Returns:
    str: The color (as a string) associated with the category that has the majority of processes associated with the gene.
    """

    # Initialize a dict to count category occurrences
    category_counts = {}
    
    # Loop over categories and processes
    for category, processes in process_dict.items():
        # Check if any process is in the gene's processes
        for process in flatten_processes(category, processes):
            if process in gene_dict.get(gene, [])['processes']:
                # If process is found, increment the count for the category
                category_counts[category] = category_counts.get(category, 0) + 1

    # Get the category with maximum count
    max_category = max(category_counts, key=category_counts.get, default=None)

    # If no category is found (i.e., the gene has no processes), return 'black'
    if max_category is None:
        return 'black'
    else:
        # Return the color for the category with maximum count
        return category_colors.get(max_category, 'black')


def flatten_processes(category, process_list):
    flat_list = []
    for process in process_list:
        if isinstance(process, dict):
            flat_list.extend(flatten_processes(category, list(process.values())))
        elif isinstance(process, list):
            flat_list.extend(flatten_processes(category, process))
        else:
            flat_list.append(process)
            
    flat_list.append(category)
    flat_list = list(set(flat_list))
    
    return flat_list