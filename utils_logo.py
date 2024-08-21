from art import *
from art import art_list, text2art


def get_colored_art(text, color_code):
    return "\n".join([color_code + line + "\033[0m" for line in text2art(text, font="small").split("\n")])


def print_side_by_side(art1, art2):
    lines1 = art1.strip().split("\n")
    lines2 = art2.strip().split("\n")

    # Pad the shorter art with empty lines if the two arts have different numbers of lines
    while len(lines1) < len(lines2):
        lines1.append("")
    while len(lines2) < len(lines1):
        lines2.append("")

    # Join the corresponding lines from the two arts with a space in between
    combined_lines = [line1 + " " + line2 for line1, line2 in zip(lines1, lines2)]

    # Print the combined lines
    print("\n".join(combined_lines))


def print_logo(msg: str):
    # ANSI color codes
    red = "\033[91m"
    magenta = "\033[95m"
    yellow = "\033[93m"
    green = "\033[92m"
    cyan = "\033[96m"

    # Generate colored art for each segment
    colored_parts = [
        get_colored_art("In", cyan),
        get_colored_art("Phe", magenta),
        get_colored_art("RNo", green),
        get_colored_art(" - ", red),
        get_colored_art("ChIP", yellow),
    ]

    # Combine the colored parts line by line
    final_art_lines = []
    for lines in zip(*[part.split("\n") for part in colored_parts]):
        final_art_lines.append("".join(lines))

    final_art = "\n".join(final_art_lines)
    print("\n", final_art)

    # Print running message
    print_side_by_side(text2art(msg, font="mini"), art("trolling"))
