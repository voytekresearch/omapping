"""WordCloud plots for OM."""

import random
from wordcloud import WordCloud

def create_wc(words_in):
    """Create WordCloud object.

    Parameters
    ----------
    words_in : list of tuple
        Words to plot, with their corresponding frequencies.

    Returns
    -------
    wc : WordCloud() object
        Wordcloud definition.
    """

    # Create the WordCloud object
    wc = WordCloud(background_color='white',
               width=800,
               height=400,
               prefer_horizontal=1,
               relative_scaling=0.5,
               min_font_size=25,
               max_font_size=80).generate_from_frequencies(words_in)

    # Change colour scheme to grey
    wc.recolor(color_func=_grey_color_func, random_state=3)

    return wc

############################################################################################
############################################################################################
############################################################################################

def _grey_color_func(word, font_size, position, orientation, random_state=None, **kwargs):
    """Function for custom coloring - use gray pallete.
    From here: https://amueller.github.io/word_cloud/auto_examples/a_new_hope.html
    """
    return "hsl(0, 0%%, %d%%)" % random.randint(25, 50)
