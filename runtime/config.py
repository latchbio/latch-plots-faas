from dataclasses import dataclass

from latch_config.config import read_config


@dataclass(frozen=True)
class Config:
    auto_reload: bool
    domain: str


config = read_config(Config)
