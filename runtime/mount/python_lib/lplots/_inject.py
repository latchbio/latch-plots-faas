from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...agent import Agent
    from ...kernel import Kernel


kernel: "Kernel" = None
agent: "Agent" = None
