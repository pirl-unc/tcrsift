# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Generate pronounceable mnemonic names from CDR3 sequences.

This module converts CDR3 amino acid sequences into human-readable names
that preserve sequence similarity - similar sequences produce similar names.
Common conserved prefixes (CASS, CAV) and suffixes (F) are stripped to focus
on the variable region that distinguishes clonotypes.

Key design choices:
- Inserted vowels use O, U, and diphthongs (ai, ei, oo, etc.) - not A, E, I, Y alone
- This allows distinguishing original amino acid vowels from inserted ones
- Different consonant cluster rules for word start, middle, and end positions
- Diphthong choice is based on surrounding consonants for natural-sounding names
"""

# Common CDR3 prefixes/suffixes (data-driven from typical TCR sequences)
DEFAULT_PREFIXES = ["CASS", "CAS", "CAV", "CA"]
DEFAULT_SUFFIXES = ["F"]

# Vowels that exist as amino acids - these come from the sequence
_AMINO_VOWELS = set("aeiy")

# Single vowels safe to insert (O and U are NOT amino acid codes)
_SAFE_SINGLES = ["o", "u"]

# Diphthongs/double vowels for insertion (more varied, natural-sounding)
# These are clearly distinguishable from single amino acid vowels
_DIPHTHONGS = [
    "ai",
    "ei",
    "oi",
    "au",
    "ou",  # Classic diphthongs
    "oo",
    "ee",  # Long vowels
    "ao",
    "eo",
    "io",  # O-ending
]

# All vowels for pronunciation purposes
_ALL_VOWELS = set("aeiouy")

# Consonant clusters acceptable at WORD START
_START_CLUSTERS = {
    "bl",
    "br",
    "ch",
    "cl",
    "cr",
    "dr",
    "dw",
    "fl",
    "fr",
    "gl",
    "gr",
    "gw",
    "kn",
    "kr",
    "ph",
    "pl",
    "pr",
    "ps",
    "qu",
    "sc",
    "sh",
    "sk",
    "sl",
    "sm",
    "sn",
    "sp",
    "st",
    "sw",
    "th",
    "tr",
    "tw",
    "wh",
    "wr",
}

# Consonant clusters acceptable in WORD MIDDLE
_MID_CLUSTERS = _START_CLUSTERS | {
    # Nasals + stops
    "mb",
    "mp",
    "nc",
    "nd",
    "ng",
    "nk",
    "nt",
    "nch",
    "nz",
    # Liquids + stops
    "lb",
    "lc",
    "ld",
    "lf",
    "lg",
    "lk",
    "lm",
    "ln",
    "lp",
    "ls",
    "lt",
    "lv",
    "rb",
    "rc",
    "rd",
    "rf",
    "rg",
    "rk",
    "rl",
    "rm",
    "rn",
    "rp",
    "rs",
    "rt",
    "rv",
    "rz",
    # Other common mid-word clusters
    "ct",
    "ft",
    "pt",
    "xt",
    "gn",
    "ks",
    "ps",
}

# Consonant clusters acceptable at WORD END
_END_CLUSTERS = {
    "ch",
    "ck",
    "ct",
    "ft",
    "fth",
    "ld",
    "lf",
    "lk",
    "lm",
    "ln",
    "lp",
    "ls",
    "lt",
    "lth",
    "lv",
    "mb",
    "mp",
    "mph",
    "nc",
    "nch",
    "nd",
    "ng",
    "ngth",
    "nk",
    "ns",
    "nt",
    "nth",
    "nz",
    "ph",
    "ps",
    "pt",
    "rb",
    "rc",
    "rd",
    "rf",
    "rg",
    "rk",
    "rl",
    "rm",
    "rn",
    "rp",
    "rs",
    "rt",
    "rth",
    "rv",
    "rz",
    "sh",
    "sk",
    "sp",
    "st",
    "sth",
    "th",
    "ts",
    "xt",
}

# Consonant-to-vowel affinity: which vowels sound good after which consonants
_CONSONANT_VOWEL_AFFINITY = {
    # Labials prefer rounded vowels
    "b": ["o", "oo", "au", "ou"],
    "p": ["o", "oo", "au", "ei"],
    "m": ["o", "oo", "ai", "au"],
    "f": ["ai", "ei", "oo", "o"],
    "v": ["ai", "ei", "o", "u"],
    "w": ["ai", "ei", "oo", "o"],
    # Dentals/alveolars
    "t": ["ai", "ei", "o", "oo"],
    "d": ["ai", "ei", "o", "oo"],
    "n": ["ai", "ei", "o", "oo"],
    "s": ["ai", "ei", "o", "oo"],
    "z": ["ai", "ei", "oo", "o"],
    "l": ["ai", "ei", "oo", "o"],
    "r": ["ai", "ei", "oo", "o"],
    # Velars prefer back vowels
    "k": ["ai", "oo", "o", "au"],
    "g": ["ai", "oo", "o", "au"],
    "h": ["ai", "ei", "oo", "o"],
    # Others
    "c": ["ai", "oo", "o", "ei"],
    "j": ["ai", "ei", "oo", "o"],
    "q": ["ai", "oo", "o", "u"],
    "x": ["ai", "ei", "o", "oo"],
    "y": ["oo", "ai", "o", "ei"],
}
_DEFAULT_VOWELS = ["o", "ai", "oo", "ei"]


def _choose_vowel(prev_consonant: str, position: int) -> str:
    """Choose a vowel to insert based on preceding consonant and position.

    Args:
        prev_consonant: The consonant before the insertion point
        position: Position in the sequence (affects variety)

    Returns:
        Vowel or diphthong to insert
    """
    options = _CONSONANT_VOWEL_AFFINITY.get(prev_consonant.lower(), _DEFAULT_VOWELS)
    return options[position % len(options)]


def _make_pronounceable(seq: str) -> str:
    """Insert vowels to make a sequence pronounceable.

    Uses O, U, and diphthongs for inserted vowels to distinguish from
    original amino acid vowels (A, E, I, Y).

    Args:
        seq: Amino acid sequence (will be lowercased)

    Returns:
        Pronounceable string with vowels inserted as needed
    """
    seq = seq.lower()
    if not seq:
        return ""

    result = []
    insert_count = 0  # Tracks number of insertions for variety

    i = 0
    while i < len(seq):
        c = seq[i]

        # If current char is a vowel, just add it
        if c in _ALL_VOWELS:
            result.append(c)
            i += 1
            continue

        # Q needs a vowel after it - use 'o' if not followed by a/e/i/o
        if c == "q":
            result.append(c)
            next_char = seq[i + 1] if i + 1 < len(seq) else ""
            if next_char not in "aeio":
                result.append("o")
                insert_count += 1
            i += 1
            continue

        # Current char is a consonant
        # Count trailing consonants in result
        trail = 0
        for x in reversed(result):
            if x in _ALL_VOWELS:
                break
            trail += 1

        # Determine position context
        at_start = not any(x in _ALL_VOWELS for x in result)

        # Decide if we need to insert a vowel
        need_insert = False

        if trail == 0:
            # Just had a vowel, single consonant is fine
            pass
        elif trail == 1:
            # One consonant trailing, check if pair is valid
            pair = result[-1] + c
            clusters = _START_CLUSTERS if at_start else _MID_CLUSTERS
            if pair not in clusters:
                need_insert = True
        elif trail >= 2:
            # Two+ consonants trailing, check if adding another is valid
            triple = "".join(result[-2:]) + c
            pair = result[-1] + c
            clusters = _START_CLUSTERS if at_start else _MID_CLUSTERS
            if triple not in clusters and pair not in clusters:
                need_insert = True

        if need_insert:
            # Choose vowel based on the last consonant and position
            prev_cons = result[-1] if result else "t"
            vowel = _choose_vowel(prev_cons, insert_count)
            result.append(vowel)
            insert_count += 1

        result.append(c)
        i += 1

    return "".join(result)


def _split_into_words(name: str, target_len: int = 4, min_word: int = 3) -> list[str]:
    """Split a pronounceable string into words.

    Splits at vowel boundaries, respecting consonant cluster rules for
    word endings and beginnings.

    Args:
        name: Pronounceable string to split
        target_len: Target length for each word
        min_word: Minimum word length

    Returns:
        List of word strings
    """
    if len(name) <= target_len + min_word:
        return [name] if name else []

    words = []
    current = ""

    i = 0
    while i < len(name):
        c = name[i]
        current += c

        # Check if we should split after this character
        if len(current) >= target_len and i < len(name) - min_word:
            remaining = name[i + 1 :]

            # Good split point: after a vowel, before a consonant
            if c in _ALL_VOWELS and remaining:
                if remaining[0] not in _ALL_VOWELS:
                    # Check if remaining starts with valid cluster
                    can_split = True
                    if len(remaining) >= 2 and remaining[1] not in _ALL_VOWELS:
                        pair = remaining[0:2]
                        if pair not in _START_CLUSTERS:
                            can_split = False

                    if can_split and len(remaining) >= min_word:
                        words.append(current)
                        current = ""

        i += 1

    if current:
        if len(current) < min_word and words:
            words[-1] += current
        else:
            words.append(current)

    return words


def _capitalize_words(words: list[str]) -> str:
    """Capitalize each word and join with spaces."""
    return " ".join(w.capitalize() for w in words if w)


def tcr_name(
    seq: str,
    prefixes: list[str] | None = None,
    suffixes: list[str] | None = None,
) -> str:
    """Generate a pronounceable mnemonic name from a CDR3 sequence.

    Similar sequences will produce similar names. Common conserved
    prefixes (CASS, CAV) and suffixes (F) are stripped to focus on
    the variable region.

    Inserted vowels use O, U, and diphthongs (ai, ei, oo, etc.), since
    A, E, I, Y are amino acid codes. This allows distinguishing original
    sequence characters from inserted ones.

    Args:
        seq: CDR3 amino acid sequence (e.g., "CASSLGQAYEQYF")
        prefixes: List of prefixes to strip (default: CASS, CAS, CAV, CA)
        suffixes: List of suffixes to strip (default: F)

    Returns:
        Human-readable name

    Examples:
        >>> tcr_name("CASSLGQAYEQYF")  # LGQ needs vowel insertion
        'Laigo Qaye Qy'
        >>> tcr_name("CASSLAGAYEQYF")  # LAGA has original A
        'Lagaye Qy'
    """
    if not seq:
        return "Anon"

    # Handle multiple sequences separated by semicolon
    if ";" in seq:
        parts = seq.split(";")
        names = [tcr_name(p.strip(), prefixes, suffixes) for p in parts if p.strip()]
        return " or ".join(names) if names else "Anon"

    seq = seq.upper().strip()

    # Use default prefixes/suffixes if not specified
    if prefixes is None:
        prefixes = DEFAULT_PREFIXES
    if suffixes is None:
        suffixes = DEFAULT_SUFFIXES

    # Strip longest matching prefix
    for prefix in sorted(prefixes, key=len, reverse=True):
        if seq.startswith(prefix):
            seq = seq[len(prefix) :]
            break

    # Strip longest matching suffix
    for suffix in sorted(suffixes, key=len, reverse=True):
        if seq.endswith(suffix):
            seq = seq[: -len(suffix)]
            break

    if not seq:
        return "Anon"

    # Make pronounceable and split into words
    pronounceable = _make_pronounceable(seq)
    words = _split_into_words(pronounceable)

    if not words:
        return pronounceable.capitalize() if pronounceable else "Anon"

    return _capitalize_words(words)
