from collections import deque


class LineDeque(object):
    """ Wrapper around collections.deque object that only adds a line to the
    deque if it's not empty.

    """

    def __init__(self):
        self._deque = deque()

    def __len__(self):
        return len(self._deque)

    def append(self, item):
        if item != '':
            self._deque.append(item)

    def appendleft(self, item):
        if item != '':
            self._deque.appendleft(item)

    def pop(self):
        return self._deque.pop()

    def popleft(self):
        return self._deque.popleft()
