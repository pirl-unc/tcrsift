def tcr_name(
        seq,
        split=True,
        extra_consonants = "BJXZ",
        extra_vowels = "OU",
        valid_vowel_clusters_start = {
            "AI",
            "AU",
            "EA",
            "EI",
            "EU",
            "EW",
            "EA",
            "EI",
            "EU",
            "IA",
            "IE",
            "OI",
            "OO",
            "OY",
        },
      valid_vowel_clusters = {
            "AI",
            "AU",
            "AY",
            "EA",
            "EE",
            "EI",
            "EU",
            "EW",
            "EY",
            "EA",
            "EE",
            "EI",
            "EU",
            "EY",
            "IA",
            "IE",
            "OI",
            "OO",
            "OU",
            "OY",
        },
        valid_consonant_clusters_start = {
            "BL",
            "BR",
            "CH",
            "CR",
            "CHR",
            "CL",
            "DR",
            "FL",

            "FR",
            "GL",

            "GR",
            "KR",
            "PH",
            "PL",
            "PR",

            "SC",
            "SCH",
            "SH",
            "SCH",
            "SL",
            "SP",
            "SPR",
            "SS",
            "ST",
            "STR",
            "TH",
            "THR",
            "TR"
        },
        valid_consonant_clusters_middle={
            "BL",
            "BR",
            "CH",
            "CR",
            "CHR",
            "CL",
            "DR",
            "FL",

            "FR",
            "FT",
            "GL",
            "GN",
            "GR",
            "GS",
            "KS",
            "LL",
            "LT",
            "NV",
            "PH",
            "PL",
            "PR",
            "RG",
            "RK",
            "RR",
            "SC",
            "SCH",
            "SH",
            "SCH",
            "SL",
            "SP",
            "SPR",
            "SS",
            "ST",
            "STR",
            "TH",
            "THR",
            "TR"
            "TT",
        },
        valid_consonant_clusters_end = {
            "CK",
            "CKS",
            "FF",
            "FT"
            "FTS",


            "LD",
            "LDS",
            "LL",

            "KS",
            "NS",
            "NG",
            "NGS",
            "NGTH",
            "PS",
            "RD",
            "RDS",
            "RP",
            "RPS",
            "RT",
            "RTS",
            "RK",

            "SH",
            "SS",
            "ST",
            "TS",
            "TT",
            "TTS",

        },
        trim_prefixes={
            "CASS": "Dr.",
            "CAS": "Prof.",
            "CSA": "Gen.",
            "CAT": "Capt.",
            "CAW": "Sir",
            "CAI": "Madame",
            "CSV": "The Honourable"
        },
        trim_suffixes={
            "QYF": "MD",
            "QFF": "PhD",
            "AFF": "Esq.",
            "LFF": "Jr.",
            "YTF": "Sr.",
        }):
    if ";" in seq:
        parts = seq.split(";")
        part_names = [
            tcr_name(part, split=True)
            for part in parts
        ]
        return " or ".join(part_names)

    name_prefix = ""
    name_suffix = ""
    if split:
        for (aa_prefix, candidate_name_prefix) in trim_prefixes.items():
            if seq.startswith(aa_prefix):
                seq = seq[len(aa_prefix):]
                name_prefix = candidate_name_prefix
                break
        for (aa_suffix, candidate_name_suffix) in trim_suffixes.items():
            if seq.endswith(aa_suffix):
                seq = seq[:-len(aa_suffix)]
                name_suffix = candidate_name_suffix
                break

        half = len(seq) // 2
        first = tcr_name(seq[:half], split=False)
        second = tcr_name(seq[half:], split=False)
        name = f"{first} {second}"
        if name_prefix:
            name = f"{name_prefix} {name}"
        if name_suffix:
            name = f"{name}, {name_suffix}"
        return name
    else:
        num_added = 0
        letters = []

        vowels_in_a_row = 0
        consonants_in_a_row = 0
        valid_clusters_start = valid_consonant_clusters_start.union(valid_vowel_clusters_start)
        valid_clusters_middle = valid_consonant_clusters_middle.union(valid_vowel_clusters)

        valid_clusters_end = valid_consonant_clusters_end.union(valid_vowel_clusters)

        for i, letter in enumerate(seq):

            last_pos = i == len(seq) - 1
            curr_is_vowel = letter in "AEIOU" or (i == last_pos and letter == "Y")
            curr_is_consonant = not curr_is_vowel
            longest_stretch = max(vowels_in_a_row, consonants_in_a_row)
            if vowels_in_a_row >= 2 and curr_is_vowel:
                letters.append(extra_consonants[num_added % len(extra_consonants)])
                num_added += 1
                consonants_in_a_row = 1
                vowels_in_a_row = 0
            elif consonants_in_a_row >= 3 and curr_is_consonant:
                letters.append(
                    extra_vowels[num_added % len(extra_vowels)] if letter != "Q" else "UO")
                num_added += 1
                consonants_in_a_row = 0
                vowels_in_a_row = 1
            elif (
                    (i == 0) or
                    (curr_is_vowel and consonants_in_a_row > 0) or
                    (curr_is_consonant and vowels_in_a_row > 0) or
                    (last_pos and i > 0 and "".join(letters[-1:] + [letter]) in valid_clusters_end) or
                    (last_pos and i > 1 and "".join(letters[-2:] + [letter]) in valid_clusters_end) or
                    (i == 1 and (letters[0] + letter) in valid_clusters_start) or
                    (i == 2 and longest_stretch == 2 and "".join(letters[-2:] + [letter]) in valid_clusters_start) or
                    (i > 1 and longest_stretch == 1 and "".join(letters[-1:] + [letter]) in valid_clusters_middle) or
                    (i > 2 and longest_stretch == 2 and "".join(letters[-2:] + [letter]) in valid_clusters_middle)
            ):
                # print(i, curr_is_vowel, last_pos, consonants_in_a_row, vowels_in_a_row, letter, letters)

                pass
            elif consonants_in_a_row > 0:
                letters.append(extra_vowels[num_added % len(extra_vowels)])
                num_added += 1
                vowels_in_a_row = 1
                consonants_in_a_row = 0
            elif vowels_in_a_row > 0:
                letters.append(extra_consonants[num_added % len(extra_consonants)])
                num_added += 1
                vowels_in_a_row = 0
                consonants_in_a_row = 1

            letters.append(letter)
            if curr_is_vowel:
                consonants_in_a_row = 0
                vowels_in_a_row += 1
            else:
                vowels_in_a_row = 0
                consonants_in_a_row += 1
    return "".join(letters).capitalize()
