import subprocess  # noqa: S404
import sys
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    import anthropic  # type: ignore


def install_and_import_anthropic() -> "anthropic":  # type: ignore
    try:
        import anthropic  # type: ignore  # noqa: PLC0415
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "anthropic"])
        import anthropic  # type: ignore  # noqa: PLC0415
    return cast("anthropic", anthropic)


anthropic = install_and_import_anthropic()
