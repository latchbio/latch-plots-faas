from typing import TypedDict


class FormInputAppearance(TypedDict, total=False):
    detail: str | None
    placeholder: str | None
    help_text: str | None
    error_text: str | None
    description: str | None


class OutputAppearance(TypedDict, total=False):
    placeholder: str | None
    outline_label: str | None
