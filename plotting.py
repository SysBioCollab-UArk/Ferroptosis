import numpy as np
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties


def calc_self_distance(kde, n_samples, x_min, x_max, num_x_pts):
    x_vals = np.linspace(x_min, x_max, num_x_pts)
    dx = np.array([x_vals[i] - x_vals[i - 1] for i in range(1, len(x_vals))])
    E_Dself = 0.5 * np.sqrt(2 / n_samples / np.pi) * np.sum(np.sqrt(kde(x_vals)[1:] * dx))
    return E_Dself


def calc_hist_distance(kde, kde_ref, x_min, x_max, num_x_pts):
    x_vals = np.linspace(x_min, x_max, num_x_pts)
    dx = np.array([x_vals[i] - x_vals[i - 1] for i in range(1, len(x_vals))])
    D = 0.5 * np.sum(np.abs(kde(x_vals)[1:] - kde_ref(x_vals)[1:]) * dx)
    return D


def write_multicolor_word(fig, x, y, word, colors, fontsize=10, fontweight='normal', additional_text=None):

    if len(colors) != len(word):
        raise Exception("Length of 'colors' (%d) does not match length of 'word' (%d)" % (len(colors), len(word)))

    # Font properties
    fontprops = FontProperties(size=fontsize, weight=fontweight)

    # Measure full word width (includes kerning)
    tp_word = TextPath((0, 0), word, prop=fontprops)
    word_width = tp_word.get_extents().width / 72 / fig.get_size_inches()[0]

    # Measure individual letter widths
    letter_widths = []
    for letter in word:
        tp_letter = TextPath((0, 0), letter, prop=fontprops)
        letter_widths.append(tp_letter.get_extents().width / 72 / fig.get_size_inches()[0])  # convert pts to figure coords

    # Estimate average inter-letter spacing from kerning
    letter_spacing = 0 if len(word) <= 1 else (word_width - sum(letter_widths)) / (len(word) - 1)

    # Draw each character
    for i, (letter, color) in enumerate(zip(word, colors)):
        fig.text(x, y, letter, color=color, ha='left', fontsize=fontsize, fontweight=fontweight)
        # Move x-position forward by width of character + additional spacing
        x += letter_widths[i] + letter_spacing

    # Add additional text, if exists
    if isinstance(additional_text, str):
        fig.text(x, y, additional_text, ha='left', fontsize=fontsize, fontweight=fontweight)
    elif isinstance(additional_text, dict):
        text = additional_text.pop('text')
        ha = additional_text.pop('ha', 'left')
        fs = additional_text.pop('fontsize', fontsize)
        fw = additional_text.pop('fontweight', fontweight)
        fig.text(x, y, text, ha=ha, fontsize=fs, fontweight=fw, **additional_text)
