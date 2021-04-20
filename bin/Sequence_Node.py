class Sequence_Node:
    seq = ""
    id = ""
    count = 0
    friends = []

    def __init__(self, seq, id, count):
        self.seq = seq
        self.id = id
        self.count = count
        self.friends = set()
        self.friends.add(self)

    def add_friend(self, friend_node):
        self.friends.update(friend_node.friends)
        self.update_friends(self.friends)

    def update_friends(self, friends_set):
        if self.friends is friends_set:
            return

        old_friends = self.friends
        self.friends = friends_set

        for friend in old_friends:
            friend.update_friends(friends_set)

    def find_best_friend(self):
        best_friend = None
        max_count = 0

        for friend in self.friends:
            if friend.count > max_count:
                best_friend = friend
                max_count = friend.count

        return best_friend
