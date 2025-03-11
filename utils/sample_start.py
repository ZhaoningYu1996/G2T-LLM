
# Sample a fragment from the fragment dictionary
import json
import random

def sample_start(start_freq, sample_method='random'):
    """
    Function to sample a fragment from the fragment dictionary
    """
    # Get the total number of fragments
    total_fragments = sum(start_freq.values())
    
    if sample_method == 'random':
        # Sample a start structure randomly
        start = random.choice(list(start_freq.keys()))
    elif sample_method == 'distributed':
        # Sample a start structure based on the distribution
        random_num = random.randint(1, total_fragments)
        for start, freq in start_freq.items():
            random_num -= freq
            if random_num <= 0:
                break
    
    return start