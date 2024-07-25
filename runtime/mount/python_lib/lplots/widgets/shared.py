from typing import Literal, TypedDict


class FormInputAppearance(TypedDict, total=False):
    detail: str | None
    placeholder: str | None
    help_text: str | None
    error_text: str | None
    description: str | None
