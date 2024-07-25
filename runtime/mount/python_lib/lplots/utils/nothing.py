from enum import Enum
from typing import ClassVar, Self


class Singleton:
    _instance: ClassVar[Self | None] = None

    def __new__(cls) -> Self:
        return cls.x

    @classmethod
    @property
    def x(cls) -> Self:
        if cls._instance is None:
            cls._instance = object.__new__(cls)

        return cls._instance


class _Nothing(Singleton):
    ...


class Nothing(Enum):
    x = _Nothing()
