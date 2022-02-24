import math
import numpy as np

from sklearn.neighbors import KDTree


class DGrid:

    def __init__(self,
                 icon_width=1,
                 icon_height=1,
                 delta=None
                 ):
        self.icon_width_ = icon_width
        self.icon_height_ = icon_height
        self.delta_ = delta

        self.grid_ = []

    def _fit(self, y):
        # calculating the bounding box
        max_coordinates = np.amax(y, axis=0)
        min_coordinates = np.amin(y, axis=0)
        bounding_box_width = max_coordinates[0] - min_coordinates[0]
        bounding_box_height = max_coordinates[1] - min_coordinates[1]
        
        # if there is no delta given, we calculate the smalles possible delta as a natural number
        if self.delta_ is None:
            print("...calc delta...")
            self.delta_ = math.ceil(math.sqrt((len(y)*self.icon_width_*self.icon_height_)/(bounding_box_height*bounding_box_width)))
        
        # defining the number of rows and columns
        nr_columns = int((self.delta_ * bounding_box_width) / self.icon_width_)
        nr_rows = int((self.delta_ * bounding_box_height) / self.icon_height_)
        print(self.delta_)
        print(nr_columns, nr_rows)

        if (nr_rows * nr_columns) < len(y):
            raise Exception("There is no space to remove overlaps! Rows: {0}, columns: {1}, data size: {2}. Try "
                            "increasing delta.".
                            format(nr_rows, nr_columns, len(y)))

        # add the original points
        def to_grid_cell(id_, x_, y_):
            return {'id': id_,
                    'x': x_,
                    'y': y_,
                    'i': 0,
                    'j': 0,
                    'dummy': False}

        for i in range(len(y)):
            self.grid_.append(to_grid_cell(i, y[i][0], y[i][1]))

        # add the dummy points
        self.add_dummy_points(nr_columns, nr_rows)

        # execute
        self.grid_ = DGrid.grid_rec(self.grid_, nr_rows, nr_columns, 0, 0)
        self.grid_.sort(key=lambda v: v.get('id'))

        transformed = []
        for i in range(len(self.grid_)):
            if self.grid_[i]['dummy'] is False:
                transformed.append(np.array([self.grid_[i]['j'] * self.icon_width_,
                                    self.grid_[i]['i'] * self.icon_height_]))

        return np.array(transformed)

    def fit_transform(self, y):
        return self._fit(y)

    def fit(self, y):
        return self._fit(y)

    @staticmethod
    def split_grid(grid, cut_point, direction):
        if direction == 'x':
            grid.sort(key=lambda cel: (cel['x'], cel['y']))
        else:
            grid.sort(key=lambda cel: (cel['y'], cel['x']))

        grid0 = grid[:cut_point]
        grid1 = []
        if cut_point < len(grid):
            grid1 = grid[-(len(grid) - cut_point):]

        return grid0, grid1

    @staticmethod
    def grid_rec(grid, r, s, i, j):
        size = len(grid)

        if size > 0:
            if size == 1:
                grid[0]['i'] = i
                grid[0]['j'] = j
            else:
                if r > s:
                    half_rows = int(math.ceil(r / 2.0))
                    grid0, grid1 = DGrid.split_grid(grid, min(size, half_rows * s), 'y')
                    DGrid.grid_rec(grid0, half_rows, s, i, j)
                    DGrid.grid_rec(grid1, (r - half_rows), s, (i + half_rows), j)
                else:
                    half_columns = int(math.ceil(s / 2.0))
                    grid0, grid1 = DGrid.split_grid(grid, min(size, half_columns * r), 'x')
                    DGrid.grid_rec(grid0, r, half_columns, i, j)
                    DGrid.grid_rec(grid1, r, (s - half_columns), i, (j + half_columns))

        return grid

    def add_dummy_points(self, nr_columns, nr_rows):
        size = len(self.grid_)
        min_x = min_y = math.inf
        max_x = max_y = -math.inf

        scatterplot = []

        for i in range(size):
            x_ = self.grid_[i]['x']
            y_ = self.grid_[i]['y']

            scatterplot.append([x_, y_])

            if min_x > x_:
                min_x = x_

            if min_y > y_:
                min_y = y_

            if max_x < x_:
                max_x = x_

            if max_y < y_:
                max_y = y_

        count = [[0] * nr_columns for _ in range(nr_rows)]

        for i in range(size):
            col = int((self.grid_[i]['x'] - min_x) / self.icon_width_)
            row = int((self.grid_[i]['y'] - min_y) / self.icon_height_)
            count[row][col] = count[row][col] + 1

        icons_area = size * self.icon_width_ * self.icon_height_
        scatterplot_area = (max_x - min_x) * (max_y - min_y)
        mask_size = int(max(3, scatterplot_area / icons_area))
        mask_size = mask_size + 1 if mask_size % 2 == 0 else mask_size

        sigma = (mask_size - 1) / 6.0
        mask = DGrid.gaussian_mask(mask_size, sigma)

        # create a kd-tree to calculate the closest points
        tree = KDTree(scatterplot, leaf_size=2)

        dummy_points_candidates = []

        for row in range(nr_rows):
            y_ = row * (max_y - min_y) / (nr_rows - 1) + min_y

            for column in range(nr_columns):
                if count[row][column] == 0:
                    density = 0

                    for i in range(mask_size):
                        for j in range(mask_size):
                            r = row - (int(mask_size / 2)) + j
                            c = column - (int(mask_size / 2)) + i

                            if (0 <= r < nr_rows) and (0 <= c < nr_columns):
                                density += mask[i][j] * count[r][c]

                    x_ = column * (max_x - min_x) / (nr_columns - 1) + min_x

                    distance, index = tree.query([[x_, y_]], 1)

                    dummy_points_candidates.append([x_, y_, density, distance[0][0]])

        dummy_points_candidates.sort(key=lambda x: (x[2], x[3]))

        nr_dummy_points = min((nr_rows * nr_columns) - size, len(dummy_points_candidates))

        for i in range(nr_dummy_points):
            self.grid_.append({'id': size + i,
                               'x': dummy_points_candidates[i][0],
                               'y': dummy_points_candidates[i][1],
                               'i': 0,
                               'j': 0,
                               'dummy': True})

        return

    @staticmethod
    def gaussian_mask(size, sigma):
        mask = [[0.0] * size for _ in range(size)]

        for i in range(size):
            y = int(i - int(size / 2))

            for j in range(size):
                x = int(j - int(size / 2))
                mask[i][j] = 1.0 / (2 * math.pi * sigma * sigma) * math.exp(-(x * x + y * y) / (2 * sigma * sigma))

        return mask
