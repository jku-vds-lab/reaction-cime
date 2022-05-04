import time

class CustomCallback():
    def __init__(self):
        pass

    def __call__(self, msg):
        pass

class TimeCallback(CustomCallback):
    def __init__(self):
        self.time = time.time()

    def __call__(self, msg):
        current_time = time.time()
        delta_time = current_time-self.time
        self.time = current_time
        print(msg, "|", "%i min %f s"%(delta_time/60, delta_time%60))