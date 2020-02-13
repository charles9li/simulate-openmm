from .line_deque import LineDeque


class NestedParser(object):
    """ This class takes in a LineDeque instance and creates a new LineDeque
    instance that is nested to reflect the nesting based on the specified
    bounding characters. The input LineDeque instance should be already created
    from the input file and cleaned of all comments.

    Parameters
    ----------
    line_deque : LineDeque
        LineDeque instance without comments
    begin_char : str
        bounding character
    end_char : str
        bounding character

    """

    # =========================================================================

    def __init__(self, line_deque, begin_char, end_char):
        self.begin_char = begin_char
        self.end_char = end_char
        self.parsed_lines = self._nested_parser(line_deque)

    # =========================================================================

    def _nested_parser(self, line_deque):
        parsed_lines, current_level = self._nested_parser_helper(line_deque, 0)
        if current_level != 0:
            raise ValueError("Enclosing character mismatch.")
        return parsed_lines

    def _nested_parser_helper(self, line_deque, current_level):
        parsed_lines = LineDeque()
        while len(line_deque) > 0:
            line = line_deque.popleft()
            bounding_char = self._determine_bounding_character(line)
            if bounding_char is not None:
                line = line.split(bounding_char, 1)
                parsed_lines.append(line[0].strip())
                line_deque.appendleft(line[1].strip())
                level_increment = self._determine_level_increment(bounding_char)
                current_level += level_increment
                if level_increment == -1:
                    inner_parsed_lines, current_level = self._nested_parser_helper(line_deque, current_level)
                    parsed_lines.append(inner_parsed_lines)
                elif level_increment == 1:
                    break
            else:
                parsed_lines.append(line)

        return parsed_lines, current_level

    # =========================================================================

    # Private helper methods for creating nested LineDeque

    def _determine_bounding_character(self, line):
        if self._has_only_begin_char(line):
            return self.begin_char
        elif self._has_only_end_char(line):
            return self.end_char
        elif self._has_both(line):
            return self._first_bounding_char(line)
        else:
            return None

    def _has_only_begin_char(self, line):
        if self.end_char in line:
            return False
        return self.begin_char in line

    def _has_only_end_char(self, line):
        if self.begin_char in line:
            return False
        return self.end_char in line

    def _has_both(self, line):
        return self.begin_char in line and self.end_char in line

    def _first_bounding_char(self, line):
        if not self._has_both(line):
            raise ValueError("Input line doesn't have both enclosing characters.")
        begin_index = line.find(self.begin_char)
        end_index = line.find(self.end_char)
        min_index = min(begin_index, end_index)
        return line[min_index]

    def _determine_level_increment(self, bounding_char):
        if bounding_char == self.begin_char:
            return -1
        elif bounding_char == self.end_char:
            return 1
        else:
            return 0
