import os
import datetime
from logging import info
from json import load


class Properties:
    def __init__(self, path):
        """
        - store data from json file
        :param path: (str) path to *.main_ctrl.json
        """

        with open(path) as file:
            self.__dict__.update(load(file, parse_int=int, parse_float=float))

        # self.auto_eval_output = True if self.eval_output == "True" else False
        path = self.working_dir = os.path.abspath(path).strip('.main_ctrl.json')

        self.scenario = os.path.basename(path)
        self.working_dir = os.path.dirname(path) + '/'  # to remove
        self.t_start = datetime.datetime.strptime(self.t_start, '%Y-%m-%d %H:%M:%S')

        info('Read inputile {}.main_ctrl.json\n    in {}'.format(self.scenario, self.working_dir))
