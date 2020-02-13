from collections import deque


class LineDeque(object):
    """ Wrapper around collections.deque object that only adds a line to the
    deque if it's not empty.

    """

    def __init__(self):
        self.deque = deque()

    def __len__(self):
        return len(self.deque)

    def append(self, item):
        if item != '':
            self.deque.append(item)

    def appendleft(self, item):
        if item != '':
            self.deque.appendleft(item)

    def pop(self):
        return self.deque.pop()

    def popleft(self):
        return self.deque.popleft()
