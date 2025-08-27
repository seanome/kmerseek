# This document is the jupyter file that shows the results of the functions made by this file
# The python file sigseq has a function called get_overlapping_kmers which outputs a list of overlaps
# This is an example of an overlap. It includes the index, hp, sequence, hashval : (30, 'phhphhhhhhhh', 'DAGDVGAAPPGA', 14427879530792476309)
# The list orders each index from smallest to largest
# This file creates a table so that the sequences are paired to eachother
# The visualization takes in the table that displays the overlaps to create the visualization
# This jupyter notebook creates a python file using magics
# The command that creates the file has been turned into a comment so that a python file is not created everytime it is ran
import pandas as pd
import seaborn as sns
from sig2kmer import degenerate_protein_chatgpt
import sourmash
import jupyter_black
import math as math
from rich.console import Console

console = Console()

jupyter_black.load()

# This function displays the aligment table

# It takes in the list created by the function get_overlapping_kmers from the Python file SigSeq,
# the names of the sequences, and the alphabet that was used
# For V2 we may want to rename this to Target and Query


import pandas as pd


def AlignmentTable(
    aligmentdatalist1: list,
    aligmentName1: str,
    aligmentdatalist2: list,
    aligmentName2: str,
    alphabet: str,
):

    columna = [
        aligmentName1 + " Index",
        alphabet,
        aligmentName1 + " protein sequence",
        "hashval",
    ]
    columnb = [
        aligmentName2 + " Index",
        alphabet,
        aligmentName2 + " protein sequence",
        "hashval2",
    ]

    sequencea = pd.DataFrame(aligmentdatalist1, columns=columna)
    sequenceb = pd.DataFrame(aligmentdatalist2, columns=columnb).drop(
        columns=["hashval2"]
    )

    return sequencea.merge(sequenceb, on=alphabet, how="left")
# This feature uses console.print from rich.console. Every match displays a different color.
# This function is used to provide a list of colors in a way that is setup for console.print
# These colors came from matplotlib's Tab20
def colors():
    darkblue = "[bold rgb(31,119,180)]"
    lightblue = "[bold rgb(174,199,232)]"
    darkorange = "[bold rgb(255,127,14)]"
    lightorange = "[bold rgb(255,187,120)]"
    darkgreen = "[bold rgb(44,160,44)]"
    lightgreen = "[bold rgb(152,223,138)]"
    darkred = "[bold rgb(214,39,40)]"
    lightred = "[bold rgb(255,152,150)]"
    darkpurple = "[bold rgb(148,103,189)]"
    lightpurple = "[bold rgb(197,176,213)]"
    darkbrown = "[bold rgb(140,86,75)]"
    lightbrown = "[bold rgb(196,156,148)]"
    darkpink = "[bold rgb(227,119,194)]"
    lightpink = "[bold rgb(247,182,210)]"
    darkgray = "[bold rgb(127,127,127)]"
    lightgray = "[bold rgb(199,199,199)]"
    darkolive = "[bold rgb(188,189,34)]"
    lightolive = "[bold rgb(219,219,141)]"
    darkturquoise = "[bold rgb(23,190,207)]"
    lightturquoise = "[bold rgb(158,218,229)]"
    color_list = [
        darkblue,
        darkorange,
        darkgreen,
        darkred,
        darkpurple,
        darkbrown,
        darkpink,
        darkgray,
        darkolive,
        darkturquoise,
        lightblue,
        lightorange,
        lightgreen,
        lightred,
        lightpurple,
        lightbrown,
        lightpink,
        lightgray,
        lightolive,
        lightturquoise,
    ]

    return color_list
# This function adds empty characters to the end of the sequence
# This allows for the sequence to be printed out in chunks of 10
# For V2 we may want to reconsider how the last 50 characters of a sequence is displayed
def corrected_length(sequence, iterations):
    length_of_sequence = len(sequence)
    total_page = 50 * iterations
    blank_letters = total_page - length_of_sequence
    sequence = sequence + " " * blank_letters
    return sequence
# This function makes the counter display for the feature 10,20,30,40,50
# Iterations means how many lines will the visualization be
# 1 iterations is 50 characters
# In V2 we may want to fix what happens at character 100


def sequence_counter(iteration):
    if iteration == 0 or iteration == 1:
        blank = " " * 8
    elif 1 < iteration < 20:
        blank = " " * 7
    elif iteration >= 20:
        blank = " " * 6
    block = iteration * 50
    return str(
        blank
        + str(10 + block)
        + " "
        + blank
        + str(20 + block)
        + " "
        + blank
        + str(30 + block)
        + " "
        + blank
        + str(40 + block)
        + " "
        + blank
        + str(50 + block)
        + " "
    )
# this function is for both the sequence and its HP
# This prints the sequences out in chunks of 10 until it hits 50 characters
def display_sequence(sequence, iteration):
    block = iteration * 50
    return (
        sequence[0 + block : 10 + block]
        + " "
        + sequence[10 + block : 20 + block]
        + " "
        + sequence[20 + block : 30 + block]
        + " "
        + sequence[30 + block : 40 + block]
        + " "
        + sequence[40 + block : 50 + block]
        + " "
    )
# The row cuts off at 50 characters
# This functions defines the cutoff points for each line
def cutoffpoints(iterations):
    cutoffs = []
    for i in range(iterations):
        cutoffs.append(51 + (i * 50))
    return cutoffs
# This created the (Position: x:y, length) part of the string that is attached to the match of the string
# This function takes in the table created by
# the table have these values matched
def end_of_match(df, row):
    column_names = list(df.columns)
    index = df.iloc[row][column_names[0]]
    sequence = df.iloc[row][column_names[2]]
    end = (
        " (Position "
        + str(index)
        + ":"
        + str(index + len(sequence))
        + ", len:"
        + str(len(sequence))
        + ")"
    )
    return end
# This splits the string if it will span over 2 rows.
# This takes in the total number of rows in the table
# This requires the uses of lists and tuples because each match needs to be the same color
# Not all matches have to be split so this keeps everything uniform and able to be tagged with the matching color


def split_match(df, iterations, row):
    column_names = list(df.columns)
    # start_position is the number that is displayed in the column __Index. For example the match could start at position 48
    start_position = df.iloc[row][column_names[4]]
    # Each row only allows for 50 characters.
    # The end position is where it starts(48) + the length of the match.
    # If the match will go beyond the allowed 50 characters it will need to be split
    end_position = start_position + len(df.iloc[row][column_names[2]])
    sequence = df.iloc[row][column_names[2]]
    splits = []
    needs_to_split = []
    for i in cutoffpoints(iterations):
        # If a cutoff point [51,101,151] is in the middle of a match it will need to be broken up
        if start_position < i < end_position:
            # This determines how many characters will be displayed in the first row
            how_many_first_row = i - start_position
            # this adds the information at the end of the match
            end_string = end_of_match(df, row)
            first_row = [start_position, (sequence[0:how_many_first_row], "")]
            second_row = [
                start_position + how_many_first_row,
                (sequence[how_many_first_row:], end_string),
            ]
            splits.append(first_row)
            splits.append(second_row)
            needs_to_split = "True"
            break
    if needs_to_split != "True":
        splits.append([start_position, (sequence, end_of_match(df, row))])
    return splits
# The sequence is displayed in chunks of 10. This function makes sure that the match has the same spacing to stay aligned
def spaces(string1, string2):
    position = 0
    string = ""
    # There is 55 characters 50+5 spaces per row
    for i in range(60):
        try:
            # If the sequence has a space, add a space
            if string1[i] == " ":
                string += " "
            else:
                string += string2[position]
                position = position + 1
                # once the match runs out of characters it means that the spaces have been fully added
        except IndexError:
            break

    return string
# This function pairs the colors properly, even if they go into differnt lines
# This also flattens out the list of lists a little bit
# In V2 we may want to incorporate this function into split_match and have the colors be added there


def colorpairs(df, iteration):
    matches = []
    for i in range(len(df)):
        matches.append(split_match(df, iteration, i))
    color_matched_matches = []

    possible_colors = colors()
    for i in range(len(matches)):
        if len(matches[i]) == 2:
            matches[i][0].append(possible_colors[i])
            matches[i][1].append(possible_colors[i])
            color_matched_matches.append(matches[i][0])
            color_matched_matches.append(matches[i][1])

        else:
            matches[i][0].append(possible_colors[i])
            color_matched_matches.append(matches[i][0])
    return color_matched_matches
# This function adds blank spaces to the matches so that they stay in alignment with the main sequence


def match_placements(df, total_iterations, current_iteration):
    color_pairs = colorpairs(df, total_iterations)
    section = []
    start_position = 0 + (current_iteration * 50)
    end_position = 50 + (current_iteration * 50)
    for i in range(len(color_pairs)):
        # If there is a match within the 50 characters of each row spaces are added to the front of the matches
        if start_position < color_pairs[i][0] < end_position:
            color_matched = [
                " " * (color_pairs[i][0] - start_position - 1) + color_pairs[i][1][0],
                color_pairs[i][1][1],
            ]
            section.append([color_matched, color_pairs[i][2]])
        else:
            continue
        # This there is not a match on this row it skips over it
    if section == []:
        section.append(["  ", "  "])
    return section
# This function is to group the matches by their colors and text
def groups(df, total_iterations, current_iteration):
    group = []
    for i in range(len(match_placements(df, total_iterations, current_iteration))):
        group.append(
            match_placements(df, total_iterations, current_iteration)[i][1]
            + spaces(
                display_sequence(ced9_seq, current_iteration),
                match_placements(df, total_iterations, current_iteration)[i][0][0],
            )
            + match_placements(df, total_iterations, current_iteration)[i][0][1]
        )
    return group
console = Console()


# mode means light or dark mode. If it is light mode text is black. If your monitor is in dark mode text is white


def kmerseekvisualization(
    aligmentdatalist1: list,
    aligmentName1: str,
    aligmentdatalist2: list,
    aligmentName2: str,
    alphabet: str,
    full_sequence: str,
    full_sequence_converted: str,
    mode: str,
):
    df = AlignmentTable(
        aligmentdatalist1,
        aligmentName1,
        aligmentdatalist2,
        aligmentName2,
        alphabet,
    )
    total_iterations = math.ceil(len(full_sequence) / 50)
    sequence = corrected_length(full_sequence, total_iterations)
    converted_sequence = corrected_length(full_sequence_converted, total_iterations)

    if mode == "dark":
        color = "[white]"

    elif mode == "light":
        color = "[black]"

    visualization_string = []
    for i in range(total_iterations):
        visualization_string.append(color + aligmentName2)
        visualization_string.append(color + sequence_counter(i))
        visualization_string.append(display_sequence(sequence, i))
        visualization_string.append(display_sequence(converted_sequence, i))
        for j in groups(df, total_iterations, i):
            visualization_string.append(j)
        visualization_string.append(" ")

    # This returns a list
    # If you want to call it and not print after, you can change it to console.print("\n".join(visualization_string))
    return visualization_string
