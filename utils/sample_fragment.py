# Sample a fragment from a fragment dictionary

import random

def sample_fragment(fragment_freq, sample_method='random'):
    """
    Function to sample a fragment from the fragment dictionary
    """
    # Get the total number of fragments
    total_fragments = sum(fragment_freq.values())
    
    if sample_method == 'random':
        fragment = random.choice(list(fragment_freq.keys()))
    elif sample_method == 'distributed':
        random_num = random.randint(1, total_fragments)
        for fragment, freq in fragment_freq.items():
            random_num -= freq
            if random_num <= 0:
                break
    
    return fragment