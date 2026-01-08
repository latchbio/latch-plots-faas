import json
from collections.abc import Mapping
from pathlib import Path
from playwright.async_api import Browser, ConsoleMessage, Page, Playwright, Request, WebSocket, async_playwright


class HeadlessBrowser:
    def __init__(self) -> None:
        self.playwright: Playwright | None = None
        self.browser: Browser | None = None
        self.page: Page | None = None

    async def start(
        self,
        notebook_url: str,
        local_storage: Mapping[str, str],
        *,
        timeout_ms: int = 30000,
    ) -> None:
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(headless=True)
        self.page = await self.browser.new_page(viewport={"width": 1280, "height": 800})

        def log_console(msg: ConsoleMessage) -> None:
            # msg.text is a property, not a method in newer Playwright versions
            print(f"[headless_console] {msg.type}: {msg.text}")

        def log_page_error(err: Exception) -> None:
            print(f"[headless_page_error] {err}")

        def log_request_failed(req: Request) -> None:
            print(f"[headless_request_failed] {req.url} {req.failure}")

        def log_ws(ws: WebSocket) -> None:
            print(f"[headless_ws_open] {ws.url}")
            ws.on("close", lambda: print(f"[headless_ws_close] {ws.url}"))
            ws.on("framereceived", lambda _: None)

        self.page.on("console", log_console)
        self.page.on("pageerror", log_page_error)
        self.page.on("requestfailed", log_request_failed)
        self.page.on("websocket", log_ws)

        storage = dict(local_storage)
        serialized = json.dumps(storage)
        await self.page.add_init_script(
            f"""
            const entries = JSON.parse({json.dumps(serialized)});
            for (const [k, v] of Object.entries(entries)) {{
                try {{
                localStorage.setItem(k, v);
                }} catch (err) {{
                console.warn("Failed to set localStorage", k, err);
                }}
            }}
            """
        )

        await self.page.goto(notebook_url, wait_until="networkidle")
        await self.page.wait_for_selector("[data-plot-ready='true']", timeout=timeout_ms)

    async def screenshot(self, path: str, *, full_page: bool = True) -> None:
        if self.page is None:
            raise RuntimeError("Headless browser page not initialized")

        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        await self.page.screenshot(path=str(p), full_page=full_page)

    async def stop(self) -> None:
        if self.browser is not None:
            await self.browser.close()
            self.browser = None
        if self.playwright is not None:
            await self.playwright.stop()
            self.playwright = None
        self.page = None
