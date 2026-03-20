from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...new_agent import Agent
    from ...kernel import Kernel


kernel: "Kernel" = None
agent: "Agent" = None
