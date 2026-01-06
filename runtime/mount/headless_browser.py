from playwright.async_api import Browser, Page, Playwright, async_playwright


class HeadlessBrowser:
    def __init__(self) -> None:
        self.playwright: Playwright | None = None
        self.browser: Browser | None = None
        self.page: Page | None = None

    async def start(self, notebook_url: str, timeout_ms: int = 30000) -> None:
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(headless=True)
        self.page = await self.browser.new_page(viewport={"width": 1280, "height": 800})

        await self.page.goto(notebook_url, wait_until="networkidle")
        await self.page.wait_for_selector("[data-agent-ready='true']", timeout=timeout_ms)

    async def stop(self) -> None:
        if self.browser is not None:
            await self.browser.close()
            self.browser = None
        if self.playwright is not None:
            await self.playwright.stop()
            self.playwright = None
        self.page = None
