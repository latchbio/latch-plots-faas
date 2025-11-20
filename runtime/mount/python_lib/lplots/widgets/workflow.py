from dataclasses import dataclass
from typing import Any, Literal, NotRequired

from latch_cli.services.launch.launch_v2 import Execution, launch

from . import _emit, _state, widget
from .button import ButtonWidget, w_button

workflow_widget_type: Literal["workflow"] = "workflow"


class WorkflowWidgetState(_emit.WidgetState[workflow_widget_type, str]):
    label: str
    readonly: bool
    wf_name: str
    params: dict[str, Any]
    version: str | None
    automatic: bool
    execution: NotRequired[Execution]


@dataclass(frozen=True, kw_only=True)
class WorkflowWidget(widget.BaseWidget):
    _button: ButtonWidget | None
    _state: WorkflowWidgetState

    @property
    def value(self) -> Execution | None:
        if not self._state.get("automatic") and self._button is not None and self._button.value:
            self._state["execution"] = launch(
                wf_name=self._state.get("wf_name"),
                params=self._state.get("params"),
                version=self._state.get("version"),
            )
            _state.submit_widget_state()

        return self._state.get("execution")

    def sample(self) -> Execution | None:
        return self._state.get("execution")


_emit.widget_registry[workflow_widget_type] = WorkflowWidget


def w_workflow(
    *,
    key: str | None = None,
    label: str,
    wf_name: str,
    version: str | None = None,
    params: dict[str, Any],
    readonly: bool = False,
    automatic: bool = False,
) -> WorkflowWidget:
    key = _state.use_state_key(key=key)

    button = None if automatic else w_button(key=f"{key}-button", label=label, readonly=readonly)
    res = WorkflowWidget(
        _key=key,
        _button=button,
        _state={
            "type": workflow_widget_type,
            "label": label,
            "readonly": readonly,
            "wf_name": wf_name,
            "params": params,
            "version": version,
            "automatic": automatic,
        },
    )

    if automatic and res._state.get("execution") is None:
        res._state["execution"] = launch(
            wf_name=wf_name,
            params=params,
            version=version,
        )

    _emit.emit_widget(key, res._state)

    return res
