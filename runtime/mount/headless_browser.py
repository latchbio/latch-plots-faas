import json
from collections.abc import Mapping
from pathlib import Path

from playwright.async_api import Browser, Page, Playwright, async_playwright


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

        self.page.on("console", lambda msg: print(f"[headless-browser-console] [{msg.type}] {msg.text}"))
        self.page.on("pageerror", lambda err: print(f"[headless-browser-error] {err}"))

        storage = dict(local_storage)
        serialized = json.dumps(storage)
        await self.page.add_init_script(
            f"""
            console.log("[HB-INIT] Starting localStorage initialization");
            const entries = JSON.parse({json.dumps(serialized)});
            console.log("[HB-INIT] Keys to set:", Object.keys(entries));
            for (const [k, v] of Object.entries(entries)) {{
                try {{
                    localStorage.setItem(k, v);
                    console.log("[HB-INIT] Set localStorage key:", k, "value length:", v.length);
                }} catch (err) {{
                    console.error("[HB-INIT] Failed to set localStorage", k, err);
                }}
            }}
            // Verify what was set
            console.log("[HB-INIT] Verifying localStorage:");
            console.log("[HB-INIT] plots.is_agent_controlled =", localStorage.getItem("plots.is_agent_controlled"));
            console.log("[HB-INIT] viewAccountId =", localStorage.getItem("viewAccountId"));
            const authData = localStorage.getItem("latch.authData");
            if (authData) {{
                try {{
                    const parsed = JSON.parse(authData);
                    console.log("[HB-INIT] latch.authData status:", parsed.status);
                    console.log("[HB-INIT] latch.authData hasAuth0Data:", !!parsed.auth0Data);
                    console.log("[HB-INIT] latch.authData hasIdToken:", !!parsed.auth0Data?.idToken);
                    console.log("[HB-INIT] latch.authData idToken prefix:", parsed.auth0Data?.idToken?.substring(0, 20));
                }} catch (e) {{
                    console.error("[HB-INIT] Failed to parse latch.authData:", e);
                }}
            }} else {{
                console.error("[HB-INIT] latch.authData not set!");
            }}
            console.log("[HB-INIT] localStorage initialization complete");
            """
        )

        await self.page.goto(notebook_url, wait_until="networkidle")
        try:
            await self.page.wait_for_selector("[data-plot-ready='true']", timeout=timeout_ms)
        except Exception:
            await self.screenshot("/var/log/headless_browser_no_selector.png")
            raise

        await self.screenshot("/var/log/headless_browser_ready.png")

    async def screenshot(self, path: str) -> None:
        if self.page is None:
            raise RuntimeError("Headless browser page not initialized")

        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        await self.page.screenshot(path=str(p), full_page=True)

    async def stop(self) -> None:
        if self.browser is not None:
            await self.browser.close()
            self.browser = None
        if self.playwright is not None:
            await self.playwright.stop()
            self.playwright = None
        self.page = None
