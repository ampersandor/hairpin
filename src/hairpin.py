"""The script contains HairpinCalculator class that calculates hairpin dg value

    .. figure:: hairpin.jpg
        :alt: hairpin dg
        :align: left
        :width: 40%

        Visualized how hairpin calculator works
"""

from datetime import datetime
import numpy as np
import re
import json
import os
from typing import Tuple

from common import sequence

class HairpinCalculator:
    """The class calculates a hairpin dg value for a sequence
    
    Args:
        minimum_bind_len (int): minimum continuous match length to form hairpin
        min_stem (int): minimum stem length
        min_loop (int): minimum loop length
        max_loop (int): maximum loop length

    Notes:
        min loop represents the minimum loop size as below::

                   .---.
                ATGC    \\
                    A    min_loop=3
            ATGCATGC    /
                   \___/
    """

    def __init__(self, minimum_bind_len: int, min_stem: int = 3, min_loop: int = 3, max_loop: int = 30):
        hairpin_json_file = os.path.join(os.path.dirname(__file__), "hairpin.json")
        with open(hairpin_json_file, "r") as json_file:
            data = json.load(json_file)
            self._stem_dg_dict = data["stem_dg_dict"]
            self._loop_dg_dict = data["loop_dg_dict"]
        self.minimum_bind_len = minimum_bind_len
        self.min_stem = min_stem
        self.min_loop = min_loop
        self.max_loop = max_loop
        self._pass = self.min_stem + (self.min_loop + 1) // 2
        self._link_test = re.compile("[|]{" + str(self.minimum_bind_len) + ",}")
        self._min_loop_wing = (min_loop + 1) // 2

    def get_stem_dg(self, seq: str):
        """It sums up all stem dg values achieved from string slicing"""
        stem_dg = 0.0
        for i in range(len(seq) - 1):
            stem_dg += self._stem_dg_dict.get(seq[i: i + 2], 0)
        return stem_dg

    def get_max_stem_dg_from_match(self, match_seq: str) -> float:
        """With the given matched nucleotides, decompose DEGs and calculate the maximum stem dG
        
        Args:
            match_seq: Matched sequences

        Returns:
            Maximum stem dG
        """
        decomposed_cal_seqs = sequence.decompose_deg(match_seq, decompose_inosine=True, decompose_n=True)
        max_stem_dg = 0
        for decomposed_cal_seq in decomposed_cal_seqs:
            max_stem_dg = max(self.get_stem_dg(decomposed_cal_seq), max_stem_dg)
        return max_stem_dg
    
    def get_max_stem_dg(self, seq_up: str, seq_dw: str):
        """
        With the given two sequences, it calculates the maximum stem dg value.
        "N" and "I" are regarded as a DEG.
        
        Args:
            seq_up: sequence as a string
            seq_dw: sequence as a string
        """
        cal_seq = str()
        for i, (up_nuc, dw_nuc) in enumerate(zip(seq_up, seq_dw)):
            if sequence.is_nucleotide(up_nuc):
                cal_seq += up_nuc
            elif sequence.is_nucleotide(dw_nuc):
                cal_seq += sequence.get_seq_rev_comp(dw_nuc)
            else:
                up_nuc_deg = set(sequence.decompose_deg(up_nuc, decompose_inosine=True, decompose_n=True))
                dw_nuc_deg = set(sequence.decompose_deg(dw_nuc, decompose_inosine=True, decompose_n=True))
                up_dw_nuc_deg = up_nuc_deg.intersection(dw_nuc_deg)
                cal_seq += up_nuc if len(up_dw_nuc_deg) > 1 else up_dw_nuc_deg.pop() if len(up_dw_nuc_deg) == 1 else ""
        return self.get_max_stem_dg_from_match(cal_seq)

    def validate_input_sequence(self, seq_len: int):
        """
        Validate that length of user input sequence is not shorter than minimum length

        Args:
            seq_len: sequence length
        """
        min_len = self.min_stem*2 + self.min_loop
        if seq_len < min_len:
            raise AttributeError(f"Minimum seq length is {min_len}, but the given seq length is {seq_len}")

    def calculate_hairpin_dg(self, seq: str) -> Tuple[bool, float, str]:
        """
        The function calculates hairpin dg value for the given sequence.

        Args:
            seq: sequence

        Returns:
            tuple containing

            - is_hairpin (bool): true if there is a predicted hairpin structure else false
            - hairpin_dg (float): Min hairpin dg if there is a result else 999.0
            - hairpin_form (str): Detail information is joined with ":;:" It consists of loop length, loop dg, \
            max stem dg, upper sequences, binding pattern and lower sequences.
            
        """
        result = dict()
        seq_len = len(seq)
        self.validate_input_sequence(seq_len)
        # 22.5.26 - If min_loop is odd, vertex is 1. Else, vertex is 0
        vertex = 1 if self.min_loop % 2 == 1 else 0
        seq_encode = sequence.encode(seq)
        seq_rev = sequence.get_seq_rev(seq)
        seq_rev_comp = sequence.encode(seq_rev, complement=True)

        # TODO: can be reduced further if we consider min_loop
        npa = np.reshape(list(seq_encode[:-self._pass]), (1, seq_len - self._pass))
        npb = np.reshape(list(seq_rev_comp[:-self._pass]), (seq_len - self._pass, 1))

        compared = np.bitwise_and(npa, npb)
        for d in range(-seq_len + self.min_stem * 2 + self.min_loop, seq_len - self.min_stem * 2 - self.min_loop + 1):
            diag = np.diagonal(compared, d)  # 짧은 stem 길이 기준으로 계산한 match 되는 nucleotide 패턴
            overlap_len = (len(diag) + self._pass) // 2  # bind length, while d representing remaining sequence length
            fold_at = overlap_len + (abs(d) + d) // 2  # add d if and only if d is positive.
            pattern = np.array(["|"] * overlap_len)  # '|' for binding representation
            total_len = overlap_len + abs(d)  # Length of longer part to check binding
            bind_pattern = "".join(np.where(diag[:overlap_len], pattern, " ")).rjust(total_len)
            # 22.5.25 - Patterns in min_loop region must be regarded as mismatch,
            #       because they could not bind because of steric hindrance
            loop_check = self._min_loop_wing - vertex
            bind_pattern = bind_pattern[:-loop_check] + " "*loop_check  # 긴 stem 기준으로 계산한 binding pattern
            seq_up = seq[:fold_at].rjust(total_len)
            seq_dw = seq_rev[: -fold_at - vertex].rjust(total_len)
            seq_vx = seq[fold_at] if vertex else "]"
            # TODO: Why calculate hairpin dG with only continuous matches? All matches should be considered
            for valid in self._link_test.finditer(bind_pattern):
                s, e = valid.span()
                # hairpin dG는 loop dG와 stem dG를 고려하여 연산
                # 따라서 match의 3' end 부분을 mismatch로 간주하여 계산해봐야 모든 케이스를 고려할 수 있음
                for temp_e in range(e, s + self.min_stem - 1, -1):
                    loop_seq_len = seq_len - 2*temp_e + abs(d)
                    if loop_seq_len > self.max_loop:
                        break
                    loop_dg = self._loop_dg_dict.get(str(loop_seq_len), 0)
                    # s, e는 긴 stem 기준으로 계산 되었으나 diag는 짧은 stem 기준이기 때문에 보정하여 match 서열 추출
                    max_stem_dg = self.get_max_stem_dg_from_match(sequence.decode(diag[s - abs(d):temp_e - abs(d)]))
                    hairpin_dg = round(loop_dg - max_stem_dg, 2)
                    result[hairpin_dg] = {
                        "seq_up": seq_up,
                        "binding": ' '*s + bind_pattern[s:temp_e] + ' '*(total_len - temp_e) + seq_vx,
                        "seq_dw": seq_dw,
                        "loop_seq_len": loop_seq_len,
                        "loop_dg": loop_dg,
                        "max_stem_dg": max_stem_dg,
                    }
            vertex = vertex ^ 1
        if not result:
            return False, 999.0, ""
        min_key = min(result)
        min_val = result[min_key]
        return (
            True,
            min_key,
            ":;:".join(
                map(
                    str,
                    [
                        min_val["loop_seq_len"],
                        min_val["loop_dg"],
                        min_val["max_stem_dg"],
                        min_val["seq_up"],
                        min_val["binding"],
                        min_val["seq_dw"],
                    ],
                )
            ),
        )


if __name__ == "__main__":
    tmp_seq = "CCCCGGAAACTGTCTGGCTGCTTTTTGCGGGG"
    start = datetime.now()
    hairpin_calculator = HairpinCalculator(minimum_bind_len=3, min_stem=3, min_loop=3)
    print(hairpin_calculator.calculate_hairpin_dg(tmp_seq))
    end = datetime.now()
    print((end - start).total_seconds())
