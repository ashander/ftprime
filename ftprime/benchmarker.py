class Timings(object):
    def __init__(self):
        self.__time_sorting = 0.0
        self.__time_appending = 0.0
        self.__time_simplifying = 0.0
        self.__time_prepping = 0.0

    @property
    def time_prepping(self):
        return self.__time_prepping

    @time_prepping.setter
    def time_prepping(self, val):
        self.__time_prepping = val

    @property
    def time_simplifying(self):
        return self.__time_simplifying

    @time_simplifying.setter
    def time_simplifying(self, val):
        self.__time_simplifying = val

    @property
    def time_sorting(self):
        return self.__time_sorting

    @time_sorting.setter
    def time_sorting(self, val):
        self.__time_sorting = val

    @property
    def time_appending(self):
        return self.__time_appending

    @time_appending.setter
    def time_appending(self, val):
        self.__time_appending = val

    @property
    def times(self):
        """
        A dict representing times spent in various
        steps.
        """
        return {'prepping': self.__time_prepping,
                'sorting': self.__time_sorting,
                'appending': self.__time_appending,
                'simplifying': self.__time_simplifying}
